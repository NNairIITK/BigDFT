!> @file
!! @author
!!    Copyright (C) 2009-2010 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
subroutine find_pfproj( Nsol,Npaw, Ngrid,Ngrid_box, rgrid, psi1s, psigrid, real_start,&
     psigrid_pseudo, dump_functions, coeffs_out  )
  use module_base
  implicit none
  !Arguments
  integer, intent(in) :: Nsol,Npaw, Ngrid, Ngrid_box,  real_start
  real(gp), intent(in) :: psi1s(Ngrid), rgrid(Ngrid)
  real(gp), intent(inout) :: psigrid(Ngrid,Nsol),psigrid_pseudo(Ngrid,Nsol)
  real(gp) , intent(out) :: coeffs_out(Npaw)
  integer :: dump_functions
  !Local variables
  real(gp) :: dumgrid(Ngrid),  dumgrid2(Ngrid), mass, mass_pseudo
  integer :: segno(Nsol), segno_pseudo(Nsol)
  integer :: i,k, igrid,j
  real(gp)  :: coeffs(Nsol), ratio

  coeffs=0.0_gp
  do i=1, Nsol-real_start+1
     do k=1, Ngrid
        dumgrid(k)=psigrid(k,i)*psi1s(k)
     enddo
     call integrate(dumgrid, dumgrid2, rgrid, Ngrid_box)
     coeffs(i)=dumgrid2(Ngrid_box)

     !! check

     dumgrid =psigrid_pseudo(:,i)*psigrid_pseudo(:,i)
     call integrate(dumgrid, dumgrid2, rgrid, Ngrid_box)
     if( abs(dumgrid2(Ngrid_box)-1.0_gp).gt.1.0D-5) Then
        print *, "  norm(psigrid_pseudo) != 1 in find_pfproj.. probably a paw dual"
        !! STOP
     endif


     dumgrid =psigrid(:,i)*psigrid(:,i)
     call integrate(dumgrid, dumgrid2, rgrid, Ngrid_box)
     if( abs(dumgrid2(Ngrid_box)-1.0_gp).gt.1.0D-5) Then
        print *, "  norm(psigrid) != 1 in find_pfproj", dumgrid2(Ngrid_box)-1.0_gp
        STOP
     endif
     !! 
     ratio = psigrid( Ngrid_box-10,i+real_start-1)/psigrid_pseudo( Ngrid_box-10,i)

     psigrid_pseudo(:,i)=ratio *psigrid_pseudo(:,i)
     
!!$
!!$     mass=0.0D0
!!$!     imass=-1
!!$     mass_pseudo=0.0D0
!!$!!!
!!$!!!     do k=Ngrid, 1,-1
!!$!!!        if( abs(psigrid(k,i))>mass) then
!!$!!!           mass= abs(psigrid(k,i))
!!$!!!           imass=k
!!$!!!        endif
!!$!!!        if ( (Ngrid-imass).gt.10 .and.  (imass-k).gt.10 ) then
!!$!!!           exit
!!$!!!        endif
!!$!!!     enddo
!!$
!!$     do k=1, Ngrid
!!$        if( abs(psigrid(k,i))>mass) mass= abs(psigrid(k,i))
!!$        if( abs(psigrid_pseudo(k,i))>mass_pseudo) mass_pseudo= abs(psigrid_pseudo(k,i))
!!$     enddo
!!$     
!!$
!!$     do k=Ngrid, 1,-1
!!$        if( abs(psigrid(k,i))>mass*0.01) then
!!$           if( psigrid(k,i).gt.0.0D0) then
!!$              segno(i)=1.0
!!$           else
!!$              segno(i)=-1
!!$           endif
!!$           exit
!!$        endif
!!$     enddo
!!$
!!$     do k=Ngrid, 1,-1
!!$        if( abs(psigrid_pseudo(k,i))>mass_pseudo*0.01) then
!!$           if( psigrid_pseudo(k,i).gt.0.0D0) then
!!$              segno_pseudo(i)=1
!!$           else
!!$              segno_pseudo(i)=-1
!!$           endif
!!$           exit
!!$        endif
!!$     end do
  enddo

  if(.true.) then
     open(unit=22,file='numerov_pfproj_ae.dat')
     do igrid=1, Ngrid
        write(22,'(200(f20.10,1x))') rgrid(igrid), (psigrid(igrid,j ), j=1,Nsol)   
     enddo
     close(unit=22)
  endif
  
  if(.true.) then
     open(unit=22,file='numerov_pfproj_ps.dat')
     do igrid=1, Ngrid
        write(22,'(200(f20.10,1x))') rgrid(igrid), (psigrid_pseudo(igrid,j ), j=1,Nsol)   
     enddo
     close(unit=22)
  endif
  
  if( Nsol .ne. Npaw .and. Npaw<Nsol-10) then
     coeffs_out=coeffs(real_start:real_start+Npaw-1)
  end if


  call  DGEMM('N','N',Ngrid ,1,   Npaw,1.0d0 ,psigrid , Ngrid ,coeffs ,Nsol, 0.0D0 , dumgrid , Ngrid)

!!$  do i=real_start,Nsol
!!$     if(dump_functions==1)  print *, coeffs(i) , coeffs(i)*segno(i)*segno_pseudo(i-real_start+1)
!!$     coeffs(i)=coeffs(i)*segno(i)*segno_pseudo(i-real_start+1)
!!$  enddo



  if(Nsol .eq. Npaw) then
     call  DGEMM('N','N',Ngrid ,1,   Nsol-real_start+1  ,1.0d0 ,psigrid_pseudo , Ngrid ,&
          coeffs(real_start) ,Nsol-real_start+1, 0.0D0 , dumgrid2 , Ngrid)
  else
     call  DGEMM('N','N',Ngrid ,1,  Npaw ,1.0d0 ,psigrid_pseudo , Ngrid ,&
          coeffs(real_start) ,Nsol-real_start+1, 0.0D0 , dumgrid2 , Ngrid)
  endif
  psigrid(:,1)=dumgrid
  psigrid(:,2)=dumgrid2

  return
END SUBROUTINE find_pfproj


subroutine find_pfproj_4tail( Nsol,Npaw, Ngrid,Ngrid_box,Ngrid_biggerbox,&
     rgrid, psi1s, psigrid, real_start,&
     ptilde, psitilde, &
     psigrid_bigger, dump_functions, coeffs_out  )
  use module_base
  implicit none
  !Arguments
  integer, intent(in) ::  Nsol,Npaw,Ngrid,Ngrid_box,Ngrid_biggerbox,real_start
  real(gp), intent(inout) :: psi1s(Ngrid), rgrid(Ngrid)
  real(gp), intent(inout) :: psigrid(Ngrid,Nsol),psigrid_bigger(Ngrid,Nsol)
  real(gp), intent(inout) :: psitilde(Ngrid,Nsol),ptilde(Ngrid,Nsol)
  real(gp) , intent(out) :: coeffs_out(Npaw)
  integer :: dump_functions
  !Local variables
  real(gp) :: dumgrid(Ngrid),  dumgrid2(Ngrid), mass, mass_pseudo
  integer :: segno(Nsol), segno_pseudo(Nsol)
  integer :: i,k, igrid,j
  real(gp)  :: coeffs(Nsol), ratio, dum, x

  print *, " in 4tail "
  !! check
  do i=1, Nsol-real_start+1

     dumgrid =psigrid_bigger(:,i)*psigrid_bigger(:,i)
     call integrate(dumgrid, dumgrid2, rgrid, Ngrid_biggerbox)
     if( abs(dumgrid2(Ngrid_biggerbox)-1.0_gp).gt.1.0D-5) Then
        print *, "  norm(psigrid_bigger) != 1 in find_pfproj_4tail"
        STOP
     endif


     dumgrid =psigrid(:,i)*psigrid(:,i)
     call integrate(dumgrid, dumgrid2, rgrid, Ngrid_box)
     if( abs(dumgrid2(Ngrid_box)-1.0_gp).gt.1.0D-5) Then
        print *, "  norm(psigrid) != 1 in find_pfproj_4tail", dumgrid2(Ngrid_box)-1.0_gp
        STOP
     endif

     dumgrid =psigrid_bigger(:,1)*psigrid(:,i)
     call integrate(dumgrid, dumgrid2, rgrid, Ngrid_box)
     print * , " >>>>>>>> " , i, "  " , dumgrid2(Ngrid_box)

  end do
  print *, " in 4tail 2,  real_start ", real_start

  do i=1, real_start-1
     print *, i
     do k=1, Ngrid
        dumgrid(k)=psigrid_bigger(k,i)*psi1s(k)
     enddo
     print *, " integ " 
     call integrate(dumgrid, dumgrid2, rgrid, Ngrid_biggerbox)
     print *, " integ OK " 

     coeffs(i)=dumgrid2(Ngrid_biggerbox)
     print *, " subtract coeff ", coeffs(i)  
     do k=1, Ngrid
        psi1s(k) = psi1s(k) -coeffs(i)*psigrid_bigger(k,i)
     enddo
  end do
  print *, " in 4tail 2.1"


  coeffs=0.0_gp
  do i=1, Nsol-real_start+1
     do k=1, Ngrid_box
        if( rgrid(k)>rgrid(Ngrid_box)*0.75_gp ) then
           x =  ( rgrid(k)-rgrid(Ngrid_box)*0.75_gp )/(  0.25_gp* rgrid(Ngrid_box) )
           dumgrid(k)=psigrid(k,i)*(psi1s(k)- psi1s(Ngrid_box) *  exp( -7.0_gp*(1.0_gp-x)**3.5_gp)   )
        else
           dumgrid(k)=psigrid(k,i)*psi1s(k)
        endif
     enddo
     call integrate(dumgrid, dumgrid2, rgrid, Ngrid_box)
     coeffs(i)=dumgrid2(Ngrid_box)
     print *, " ============= coeff ", i , " " , coeffs(i)
  end do
 
  print *, " in 4tail 3"


  print *, " Npaw ", Npaw

  coeffs_out(:)=coeffs(real_start:real_start+Npaw-1)


  do i=1, Nsol-real_start+1
     ratio = psigrid( Ngrid_box-10,i+real_start-1)/psitilde( Ngrid_box-10,i)
     print *, "psigrid ", i , "would require q correction factor ", ratio 
     !!$ psitilde(:,i)=ratio *psitilde(:,i)
  enddo
  print *, " in 4tail 4"

  dumgrid (:) = psi1s 
  dumgrid2(:) = psi1s

  do k=1, Ngrid_box
     if( rgrid(k)>rgrid(Ngrid_box)*0.75_gp ) then
        x =  ( rgrid(k)-rgrid(Ngrid_box)*0.75_gp )/(  0.25_gp* rgrid(Ngrid_box) )
        dumgrid2(k)= psi1s(Ngrid_box) *  exp( -7.0_gp*(1.0_gp-x )** 3.5_gp )   
     else
        dumgrid2(k)=0.0_gp
     end if
  enddo
     

  do i=1, Npaw
!!$     do k=1, Ngrid
!!$        dumgrid(k)=ptilde(k,i)*psi1s(k)
!!$     enddo
!!$     call integrate(dumgrid, dumgrid2, rgrid, Ngrid_box)
!!$     dum =dumgrid2(Ngrid_box)
     do k=1, Ngrid_box
        !! psi1s(k)     = 0.0_wp
        !! psi1s(k)     =   psi1s(k) -dum*psitilde(k,i)
        !! psi1s(k)     =   psi1s(k) -coeffs(i+real_start-1)*psigrid(k,i+real_start-1)
        dumgrid2(k)  =  dumgrid2(k)   +  coeffs_out(i)*psitilde(k,i)
     enddo
  end do

  print *, " in 4tail 5"

  psigrid(:,1)=dumgrid
  psigrid(:,2)=dumgrid2

  print *, " in 4tail 6"

  return
END SUBROUTINE find_pfproj_4tail


subroutine find_Scoeffs_grid( ng,  expo, Ngrid, rgrid, psi1s , gcoeffs , l )
  use module_base
  implicit none
 
  integer, intent(in) :: ng,  l, Ngrid
  real(gp), intent(in) :: rgrid(Ngrid), psi1s(Ngrid)
  real(gp), dimension(ng+1), intent(in) :: expo
  real(gp), dimension(0:ng), intent(out) :: gcoeffs

  ! -----------------------------------------------------------
  real(gp)  ::   Soverlap(0:ng,0:ng)
  real(gp)  ::   Score(0:ng), dumgrid(Ngrid), dumgrid2(Ngrid)
  integer :: i,j,k,n,INFO, LWORK
  real(gp) :: b1,b2, B, spi, pi
  real(gp) :: W(0:ng), WORK(3*(ng+1)*(ng+1))
  real(gp) ::  sum, totalpow, ggg, gamma_restricted

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
        
        ggg= gamma_restricted( (1.0_gp+totalpow)/2.0_gp   )
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
END SUBROUTINE find_Scoeffs_grid


subroutine dump_1gauwf_on_radgrid(prefix, ng , expo,psi   ,lpow   )
  use module_base
  implicit none
 
  character(*) , intent(in) ::  prefix
  integer, intent(in) :: ng,lpow
  real(gp), dimension(ng+1), intent(in) :: expo
  real(gp), dimension(0:ng), intent(in) :: psi

  ! local
  integer, parameter :: n_int=1000
  character(len=200) :: filename
  integer :: i,ig
  real(kind=8) :: r,sum

  write(filename,'(a)') prefix


  open(unit=22,file=filename)
  do i=1, 20000
     r=0.001*i
     sum=0.0
     do ig = 0,ng
        sum=sum+psi(ig)*exp( -r*r/2.0/expo(ig+1)/expo(ig+1) )
     enddo
     write(22,*) r, sum*(r**lpow)
  enddo
  close(unit=22)

END SUBROUTINE dump_1gauwf_on_radgrid


function value_at_r(r, ng , expo,psi     )
  use module_base, only: gp

  implicit none
 
  real(gp) , intent(in) ::  r
  integer, intent(in) :: ng
  real(gp), dimension(ng+1), intent(in) :: expo
  real(gp), dimension(0:ng), intent(in) :: psi

  ! local
  integer, parameter :: n_int=1000

  integer ig
  real(gp) sum
  real(gp) :: value_at_r

  sum=0.0
  do ig = 0,ng
     sum=sum+psi(ig)*exp( -r*r/2.0/expo(ig+1)/expo(ig+1) )
  enddo
  value_at_r=sum

end function value_at_r


subroutine dump_gauwf_on_radgrid(prefix, ng, noccmax, lmax, expo, psi)
  use module_base, only: gp
  implicit none

  !Arguments
  character(*) , intent(in) ::  prefix
  integer, intent(in) :: ng,noccmax, lmax
  real(gp), dimension(ng+1), intent(in) :: expo
  real(gp), dimension(0:ng,noccmax,lmax+1), intent(in) :: psi

  !Local variables
  integer, parameter :: n_int=1000
  character(len=200) :: filename
  integer :: l,i,k,ig
  real(kind=8) :: r,sum

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
END SUBROUTINE dump_gauwf_on_radgrid

subroutine dump_real_on_radgrid(prefix, ng, noccmax, lmax, expo, psi, rgrid)
  use module_base, only: gp
  implicit none

  !Arguments
  character(*) , intent(in) ::  prefix
  integer, intent(in) :: ng,noccmax, lmax
  real(gp), dimension(ng+1), intent(in) :: expo
  real(gp), dimension(0:ng,noccmax,lmax+1), intent(in) :: psi
  real(gp), dimension(ng) :: rgrid

  !Local variables
  integer, parameter :: n_int=1000
  character(len=200) :: filename
  integer :: l,i,k,ig
  real(kind=8) :: r,sum

  do i=1,noccmax
     do l=0,lmax
  
        write(filename,'(a,a1,i1,a1,i1)') prefix,'_',i,'_',l
        
        open(unit=22,file=filename)
        do k=1 ,ng
           r=rgrid(k)
           write(22,*) r, psi(k,i,l+1)
        enddo
        close(unit=22)
     enddo
  enddo
return
END SUBROUTINE dump_real_on_radgrid


subroutine abs_generator_modified(iproc,izatom,ielpsp,psppar,npspcode,ng, noccmax, lmax ,expo,&
     psi, aeval, occup, psp_modifier, &
     Nsol, Labs, Ngrid,Ngrid_box, Egrid,  rgrid , psigrid, Npaw,  PAWpatch , psipsigrid )

  use module_base, only: gp, memocc,ndebug
  use module_interfaces
  implicit none
  integer, intent(in) :: iproc,izatom,ielpsp,ng,npspcode,noccmax, lmax, Nsol, labs, Ngrid,  Ngrid_box

  real(gp), dimension(0:4,0:6), intent(in) :: psppar
  !! real(gp), pointer, intent(in) :: psppar(:,:)
  
  integer, intent(in) :: psp_modifier, Npaw
  
  real(gp), dimension(ng+1), intent(out) :: expo

  integer, parameter :: n_int=1000

  real(gp), dimension(0:ng,noccmax,lmax+1), intent(out) :: psi, Egrid(Nsol),&
       rgrid(Ngrid), psigrid(Ngrid,Nsol  )
  real(gp),   intent(out), optional  :: psipsigrid(Ngrid,Nsol  )
  real(gp), dimension(noccmax,lmax+1  ), intent(out) ::  aeval,occup
  real(gp):: PAWpatch(Npaw,Npaw)

  !local variables
  character(len=*), parameter :: subname='abs_generator_modified'
  character(len=2) :: symbol
  real(gp), parameter :: fact=4.0_gp
  real(gp), dimension(6,4) :: neleconf
  real(gp), dimension(3) :: gpot
  real(gp), dimension(6) :: ott
  real(gp), dimension(noccmax,lmax+1) ::chrg,res
  real(gp), dimension(:), allocatable :: xp,alps
  real(gp), dimension(:,:), allocatable :: vh,hsep,ofdcoef

  real(gp), dimension(:,:,:,:), allocatable :: rmt
  integer :: lpx,nsccode,mxpl,mxchg
  integer :: l,i,iocc,i_all,i_stat,  j 
  real(gp) :: alpz,alpl,rcov,rprb,zion,rij,a,a0,a0in,tt,ehomo
  real(gp) :: value_at_r,amu
  integer :: igrid, isol
  logical :: pawisactive
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
  allocate(alps(lpx+1+ndebug),stat=i_stat)
  call memocc(i_stat,alps,'alps',subname)
  allocate(hsep(6,lpx+1+ndebug),stat=i_stat)
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
     allocate(ofdcoef(3,4+ndebug),stat=i_stat)
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
     call modified_eleconf(izatom,symbol,rcov,rprb,ehomo,neleconf,nsccode,mxpl,mxchg)
  else
     call eleconf(izatom,ielpsp,symbol,rcov,rprb,ehomo,neleconf,nsccode,mxpl,mxchg,amu)
  endif

  occup(:,:)=0.000000000_gp
   do l=0,lmax
     iocc=0
     do i=1,6
        ott(i)=neleconf(i,l+1)
        if (ott(i) > 0.0_gp) then
           iocc=iocc+1
            if (iocc > noccmax) stop 'abs_generator_modified: noccmax too small'
           occup(iocc,l+1)=ott(i)
        endif
     end do

  end do



  if (psp_modifier/=1 .and. ng<1000 ) then  ! .and. psp_modifier/=-1
     
     !allocate arrays for the gatom routine
     allocate(vh(4*(ng+1)**2,4*(ng+1)**2+ndebug),stat=i_stat)
     call memocc(i_stat,vh,'vh',subname)

     allocate(xp(0:ng+ndebug),stat=i_stat)
     call memocc(i_stat,xp,'xp',subname)
     allocate(rmt(n_int,0:ng,0:ng,lmax+1+ndebug),stat=i_stat)
     call memocc(i_stat,rmt,'rmt',subname)

     !can be switched on for debugging
     !if (iproc.eq.0) write(*,'(1x,a,a7,a9,i3,i3,a9,i3,f5.2)')&
     !     'Input Guess Generation for atom',trim(atomname),&
     !     'Z,Zion=',izatom,ielpsp,'ng,rprb=',ng+1,rprb

     rij=3._gp
     ! exponents of gaussians
     ! print *, " ESPONENTI " 
     ! print *, " alpz " , alpz
     a0in=alpz
     a0=a0in/rij
     !       tt=sqrt(sqrt(2._gp))
     tt=2._gp**(.3_gp)
     do i=0,ng
        a=a0*tt**i
        xp(i)=.5_gp/a**2
        ! print *, " xp(", i,")", xp(i)
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
  else
     psi=0.0_gp
  endif

!!!  call gatom(rcov,rprb,lmax,lpx,noccmax,occup,&
!!!       zion,alpz,gpot,alpl,hsep,alps,vh,xp,rmt,fact,n_int,&
!!!       aeval,ng,psi,res,chrg)

  if(psp_modifier==-1) then
     if(iproc==0) print *, " calling gatom_modified_eqdiff, Labs",Labs

     if(ng <1000) then

        call gatom_modified_eqdiff(rcov,rprb,lmax,lpx,noccmax,occup,&
             zion,alpz,gpot,alpl,hsep,alps,vh,xp,rmt,fact,n_int,&
             aeval,ng,psi,res,chrg,&
             Nsol, Labs, Ngrid,Egrid,  rgrid , psigrid,Npaw, PAWpatch )
     else
        call atom_numeric_eqdiff(rprb,lmax,noccmax,occup,&
             zion,&
             aeval,psi,&
             Nsol, Labs, Ngrid,Egrid,  rgrid , psigrid )
     endif


     stop
  else 
     if (psp_modifier==1) then
        call atom_numeric(rprb,lmax,noccmax,occup,&
             zion,&
             aeval,psi,&
             Nsol, Labs, Ngrid,Ngrid_box, Egrid,  rgrid , psigrid )
     else
        if(.not. present(psipsigrid) ) then
           call gatom_modified(rcov,rprb,lmax,lpx,noccmax,occup,&
                zion,alpz,gpot,alpl,hsep,alps,vh,xp,rmt,fact,n_int,&
                aeval,ng,psi,res,chrg,&
                Nsol, Labs, Ngrid,Ngrid_box,Egrid,  rgrid , psigrid,Npaw,  PAWpatch )
        else
           print *, "chiamo gatom_modified con psipsigrid "
           PAWpatch=0.0_gp
           call gatom_modified(rcov,rprb,lmax,lpx,noccmax,occup,&
                zion,alpz,gpot,alpl,hsep,alps,vh,xp,rmt,fact,n_int,&
                aeval,ng,psi,res,chrg,&
                Nsol, Labs, Ngrid,Ngrid_box,Egrid,  rgrid , psigrid,Npaw,  PAWpatch,&
                psipsigrid)           
        endif
     endif
     !! the operation below ensure that AE and pseudo wf 
     !! coincides for big r
     !! In the case of paw fitting, however, psigrid is the dual
     !! ptilde, while psitilde has been fitted to psigrid from input
     !! In this latter case the operation below is no more necessary
     !! ( Moreover ptilde might have another  sign than psitilde due the duality operation ?)
     pawisactive=.false.
     do i=1, Npaw
        do j=1, Npaw
           if( PAWpatch(i,j) /= 0.0_gp ) then
              pawisactive= .true.
           endif
        end do
     end do
     if( .not. pawisactive) then
        do isol=1,nsol
           if( psigrid(Ngrid_box-1, isol)<0) then
              do igrid=1, ngrid
                 psigrid(igrid, isol)=-psigrid(igrid, isol)
              end do
           endif
        enddo
     endif
  endif

  if (psp_modifier/=1 ) then
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

     i_all=-product(shape(xp))*kind(xp)
     deallocate(xp,stat=i_stat)
     call memocc(i_stat,i_all,'xp',subname)

     i_all=-product(shape(rmt))*kind(rmt)
     deallocate(rmt,stat=i_stat)
     call memocc(i_stat,i_all,'rmt',subname)
  endif


  i_all=-product(shape(hsep))*kind(hsep)
  deallocate(hsep,stat=i_stat)
  call memocc(i_stat,i_all,'hsep',subname)
  i_all=-product(shape(alps))*kind(alps)
  deallocate(alps,stat=i_stat)
  call memocc(i_stat,i_all,'alps',subname)

END SUBROUTINE abs_generator_modified



subroutine iguess_generator_modified(izatom,ielpsp,zion,psppar,npspcode,ng,nl,&
     nmax_occ,noccmax,lmax,occup,expo,psiat,enlargerprb, gaenes_aux)
  use module_base
  implicit none
  logical, intent(in) :: enlargerprb
  integer, intent(in) :: ng,npspcode,nmax_occ,lmax,noccmax,ielpsp,izatom
  real(gp), intent(in) :: zion
  integer, dimension(lmax+1), intent(in) :: nl
  real(gp), dimension(0:4,0:6), intent(in) :: psppar
  real(gp), dimension(noccmax,lmax+1), intent(in) :: occup
  real(gp), dimension(ng+1), intent(out) :: expo
  real(gp), dimension(ng+1,nmax_occ), intent(out) :: psiat
  real(gp), dimension(nmax_occ) :: gaenes_aux


  !local variables
  character(len=*), parameter :: subname='iguess_generator'
  integer, parameter :: n_int=100
  real(gp), parameter :: fact=4.0_gp
  character(len=2) :: symbol
  integer :: lpx,nsccode,mxpl,mxchg
  integer :: l,i,j,iocc,i_all,i_stat
  real(gp) :: alpz,alpl,amu,rprb,rij,a,a0,a0in,tt,ehomo,rcov

  real(kind=8), dimension(6,4) :: neleconf
  real(gp), dimension(4) :: gpot
  real(gp), dimension(noccmax,lmax+1) :: aeval,chrg,res
  real(gp), dimension(:), allocatable :: xp,alps
  real(gp), dimension(:,:), allocatable :: vh,hsep,ofdcoef
  real(gp), dimension(:,:,:), allocatable :: psi
  real(gp), dimension(:,:,:,:), allocatable :: rmt

  !filename = 'psppar.'//trim(atomname)

  lpx=0
  lpx_determination: do i=1,4
     if (psppar(i,0) == 0.0_gp) then
     exit lpx_determination
     else
        lpx=i-1
     end if
  end do lpx_determination

  allocate(alps(lpx+1+ndebug),stat=i_stat)
  call memocc(i_stat,alps,'alps',subname)
  allocate(hsep(6,lpx+1+ndebug),stat=i_stat)
  call memocc(i_stat,hsep,'hsep',subname)

  !assignation of radii and coefficients of the local part
  alpz=psppar(0,0)
  alpl=psppar(0,0)
  alps(1:lpx+1)=psppar(1:lpx+1,0)
  gpot(1:4)=psppar(0,1:4)

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
     allocate(ofdcoef(3,4+ndebug),stat=i_stat)
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

  !!Just for extracting the covalent radius and rprb
  call eleconf(izatom,ielpsp,symbol,rcov,rprb,ehomo,neleconf,nsccode,mxpl,mxchg,amu)

  if (enlargerprb) then
     !experimental
     rprb=100.0_gp
  end if

!  occup(:,:)=0.0_gp
!   do l=0,lmax-1
!     iocc=0
!     do i=1,6
!        if (elecorbs(i,l+1) > 0.0_gp) then
!           iocc=iocc+1
!           !print *,'elecorbs',i,l,elecorbs(i,l+1),noccmax
!            if (iocc > noccmax) stop 'iguess_generator: noccmax too small'
!           occup(iocc,l+1)=elecorbs(i,l+1)
!        endif
!     end do
!     nl(l+1)=iocc
!  end do

  !allocate arrays for the gatom routine
  allocate(vh(4*(ng+1)**2,4*(ng+1)**2+ndebug),stat=i_stat)
  call memocc(i_stat,vh,'vh',subname)
  allocate(psi(0:ng,noccmax,lmax+ndebug),stat=i_stat)
  call memocc(i_stat,psi,'psi',subname)
  allocate(xp(0:ng+ndebug),stat=i_stat)
  call memocc(i_stat,xp,'xp',subname)
  allocate(rmt(n_int,0:ng,0:ng,lmax+ndebug),stat=i_stat)
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
  do l=0,lmax-1
     do iocc=1,noccmax
        do i=0,ng
           psi(i,iocc,l+1)=0.0_gp
        end do
     end do
  end do

  call crtvh(ng,lmax-1,xp,vh,rprb,fact,n_int,rmt)

  call gatom(rcov,rprb,lmax-1,lpx,noccmax,occup,&
       zion,alpz,gpot,alpl,hsep,alps,vh,xp,rmt,fact,n_int,&
       aeval,ng,psi,res,chrg)

  !post-treatment of the inguess data
  do i=1,ng+1
     expo(i)=sqrt(0.5_gp/xp(i-1))
  end do

  i=0
  do l=1,4
     do iocc=1,nl(l)
        i=i+1
        !occupat(i)=occup(iocc,l)
        do j=1,ng+1
           psiat(j,i)=psi(j-1,iocc,l)
        end do
        gaenes_aux(i) = aeval(iocc,l)
     end do
  end do

  i_all=-product(shape(vh))*kind(vh)
  deallocate(vh,stat=i_stat)
  call memocc(i_stat,i_all,'vh',subname)
  i_all=-product(shape(psi))*kind(psi)
  deallocate(psi,stat=i_stat)
  call memocc(i_stat,i_all,'psi',subname)
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

END SUBROUTINE iguess_generator_modified


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
END SUBROUTINE integrate


function pow(x,n)
  use module_base, only: gp,wp
  implicit none
  real(gp), intent(in) :: x
  integer, intent(in) :: n
  real(gp) :: pow
  pow=x**n
end function pow

function phase(E, N, rgrid, V, nonloc, y, l, normalize, onlyout)
  use module_base, only: gp,wp
  implicit none
  !Arguments
  integer :: N, normalize, onlyout
  real(gp) :: E,rgrid(N),V(N), nonloc(N), y(N),l
  real(gp) :: phase
  !Local variables
  integer :: ii, i,j
  real(gp) :: ypi
  
  integer :: yNcross, yflag
  real(gp)  :: dh,dl,dh2,dl2,dp,dl3,dhdl2,dh2dl,dh3,add1,add2,deno,num
  real(gp)  :: Gh,Gl,G0,Nlc,Nla,Nlb
  
  real(gp) :: Ga,Gb,Gc
  real(gp) :: ya,yb,yc
  real(gp) :: r, PI
  
  real(gp) :: fact,norm, func(N), funcdum(N), pow
  integer :: count

  PI=4*atan(1.0d0)
  
  phase=0
  yNcross=0
  yflag=0
  
  ! -- Calcul du point intermediaire o le V est minimum
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
        ! stop
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
        
        ! // **
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
  
  ! // **
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
  
  ! //------------ Propagation de I rinf --------------------------
  
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
        
        ! // **
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

  ! // **
  !    // ** Termine aggiuntivo
  !    // **
  !    // 6*l^3*Nlc + 
  !    // Nla*(-12*h^2*l - 6*h*l^2 + h^3*l*(h + l)*Gc) + 
  !    // Nlb*(-24*h^2*l - 30*h*l^2 - 6*l^3 +      2*h^3*l*(h + l)*Gc + 3*h^2*l^2*(h + l)*Gc)
  !    // / al denominatore
  !    //
  !    //  6*h*(h + l)*(-12 + (h^2 + h*l - l^2)*Gc)

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
  integer :: i
  real(gp) :: Phase, but
  real(gp) :: PI

  PI=4.0*atan(1.0)
  Elow=0;

  do i=1, NGRID
      if(V(i) .lt. Elow) Elow=V(i)
   enddo
  if( Elow .lt. -Z*Z) Elow = -Z*Z
  Elow=-Z*Z

  Ehigh = 0.0;

  pathh = Phase(Ehigh,NGRID,r,v,nonloc,y,  l ,0, 0);
  pathl = Phase(Elow ,NGRID, r,v,nonloc,y,  l ,0, 0);
 

!!!  print *, Ehigh, pathh
!!!  print *, Elow, pathl

  but= PI*(nsol-l-1)

!!!  print *, " Z " , Z

  if( pathl.gt.but)  then
     print *, " pathl>but " 
     print *, " Elow " , Elow
     print *, " now exiting , routine schro" 
     stop
  endif
  
      
  if(but .gt. pathh) then
     Ehigh = (but+1)*(but+1)/r(NGRID-1)/r(NGRID-1)
     do while( but .gt. Phase(Ehigh ,NGRID, r,V,nonloc,y,  l , 0, 0 ) )
         Ehigh =2*Ehigh;
      enddo
   endif

   do while(1.gt.0) 

      Eguess = (Elow+Ehigh)/2

      fase=  Phase(Eguess ,NGRID, r,V,nonloc,y,  l , 0, 0)

      if( fase.gt.but) then
         Ehigh=Eguess
      else
         Elow=Eguess
      endif

      if(dabs(but-fase).lt.1.0e-8) exit
   enddo

   
   fase  = Phase(Eguess,NGRID, r,v,nonloc,y, l,1 ,0)
   E=Eguess;

   return

END SUBROUTINE schro


subroutine gatom_modified(rcov,rprb,lmax,lpx,noccmax,occup,&
                 zion,alpz,gpot,alpl,hsep,alps,vh,xp,rmt,fact,nintp,&
                 aeval,ng,psi,res,chrg,&
                 Nsol, Labs, Ngrid,Ngrid_box, Egrid,  rgrid , psigrid, Npaw, PAWpatch, &
                 psipsigrid)
  use module_base, only: gp
  use esatto

  implicit real(gp) (a-h,o-z)
  logical :: noproj, readytoexit
  integer, parameter :: n_int=1000
  dimension psi(0:ng,noccmax,lmax+1),aeval(noccmax,lmax+1),&
       hh(0:ng,0:ng),ss(0:ng,0:ng),eval(0:ng),evec(0:ng,0:ng),&
       gpot(3),hsep(6,lpx+1),rmt(n_int,0:ng,0:ng,lmax+1),&
       pp1(0:ng,lpx+1),pp2(0:ng,lpx+1),pp3(0:ng,lpx+1),alps(lpx+1),&
       potgrd(n_int),&
       rho(0:ng,0:ng,lmax+1),rhoold(0:ng,0:ng,lmax+1),xcgrd(n_int),&
       occup(noccmax,lmax+1),chrg(noccmax,lmax+1),&
       vh(0:ng,0:ng,4,0:ng,0:ng,4),&
       res(noccmax,lmax+1),xp(0:ng),& 
       psigrid(Ngrid, Nsol),psigrid_naked(Ngrid,Nsol),&
       psigrid_naked_2(Ngrid,Nsol), projgrid(Ngrid,3), &
       rhogrid(Ngrid), potgrid(Ngrid), psigrid_not_fitted(Ngrid,Nsol),&
       psigrid_not_fitted_2(Ngrid,Nsol),&
       vxcgrid(Ngrid), &
       Egrid(nsol), ppgrid(Nsol,3), work(nsol*nsol*2), &
       H(Nsol, Nsol), &
       H_2(Nsol, Nsol), &
       Hcorrected(Nsol, Nsol), &
       Hadd(Nsol, Nsol), Egrid_tmp(Nsol),Egrid_tmp_2(Nsol), Etofit(Nsol), &
       Soverlap(Nsol,Nsol), Tpsigrid(Nsol,Ngrid ),Tpsigrid_dum(Nsol, Ngrid),valuesatp(Nsol), &
       PAWpatch(Npaw, Npaw ), Spsitildes(Npaw, Npaw), genS(Nsol,Nsol), genH(Nsol,Nsol) , dumH(Nsol,Nsol)

  real(gp) , optional :: psipsigrid(Ngrid, Nsol)
  

  real(gp) :: rgrid(Ngrid), ene_m, ene_p, factadd, rcond, fixfact
  real(gp), target :: dumgrid1(Ngrid),dumgrid2(Ngrid), dumgrid3(Ngrid)
  logical dofit
  integer real_start, iocc, iwork(Nsol), INFO, volta, ngrid_box_2
  character(1) EQUED
  integer ipiv(Nsol), Npaw


  dofit=.false.
  if( Egrid(1).ne.0.0_gp .and. Egrid(2) .ne.0.0_gp  .and. Npaw.ne.Nsol) then
     dofit=.true.
     Etofit=Egrid
     do isol=1, Nsol
        valuesatp(isol)= psigrid(Ngrid_box-10,isol)
     enddo
  endif

  if (nintp.ne.n_int) stop 'n_int><nintp xabs'

  do l=0,lmax
     if (occup(1,l+1).gt.0._gp) lcx=l
  end do
  !write(6,*) 'lcx',lcx
 
  noproj=.true.
  do l=1,lpx+1
     noproj = noproj .and. (alps(l) .eq. 0._gp)
  end do


  

! projectors, just in case
  if ( .not. noproj) then
     do l=0,lpx
        gml1=sqrt( gamma_restricted(real(l,gp)+1.5_gp) / (2._gp*alps(l+1)**(2*l+3)) )
        gml2=sqrt( gamma_restricted(real(l,gp)+3.5_gp) / (2._gp*alps(l+1)**(2*l+7)) )&
            /(real(l,gp)+2.5_gp)
        gml3=sqrt( gamma_restricted(real(l,gp)+5.5_gp) / (2._gp*alps(l+1)**(2*l+11)) )&
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
                 projgrid(igrid,iorder) = (  exp( -r*r*tt)*(r**l *r) )* r**(2*(iorder-1))
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
  readytoexit=.false.
  big_loop: do it=1,50
     if(it.eq.50) readytoexit=.true.
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


     if( readytoexit) then

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
     endif
  
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


     if(readytoexit) then
        
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
     endif
! ------------------------------------------------------------------------------------



 

     loop_l: do l=0,lmax
        gml=.5_gp*gamma_restricted(.5_gp+real(l,gp))

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
              hh(i,j)=hh(i,j)+ gpot(1)*.5_gp*gamma_restricted(1.5_gp+real(l,gp))*tt**(1.5_gp+real(l,gp))&
                   + (gpot(2)/alpl**2)*.5_gp*gamma_restricted(2.5_gp+real(l,gp))*tt**(2.5_gp+real(l,gp))&
                   + (gpot(3)/alpl**4)*.5_gp*gamma_restricted(3.5_gp+real(l,gp))*tt**(3.5_gp+real(l,gp))
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
     ! write(6,*) 'evdiff',it,tt, readytoexit
     if (tt.lt.1.e-12_gp) then
        if( readytoexit) then
           exit big_loop
        endif
        readytoexit=.true.
     end if
  end do big_loop
! End of the big loop


  do igrid=1, ngrid
     r=rgrid(igrid)
     potgrid(igrid)=potgrid(igrid)+ 0.5_gp*labs*(labs+1.0_gp)/r/r
  enddo
  
  dumgrid1(:)=0.0_gp
  do isol=1,nsol
      psigrid_naked(:,isol)=0.0_gp
     call schro(Egrid(isol),rgrid,potgrid,dumgrid1,psigrid_naked(1,isol),ngrid_box,isol+labs,labs*1.0_gp,zion)
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

  Rbox=rgrid(Ngrid_box)
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

        do igrid=1,Ngrid_box
           r=rgrid(igrid)
           rr = r/alpz
           rrb=r/Rbox
           dumgrid1(igrid)=psigrid_naked(igrid,i)*psigrid_naked(igrid,j) *exp(- 0.5_gp * rr**2 )
           !! 
           !!  *(  1.0-2*rrb+rrb**2 ) 
        enddo

        call integrate(dumgrid1,dumgrid2,rgrid,Ngrid_box)
        Hadd(i,j)=dumgrid2(Ngrid_box)
     enddo
  enddo

!!$  if(present(psipsigrid)) then
!!$     ngrid_box_2=psipsigrid(1,1)
!!$     dumgrid1(:)=0.0_gp
!!$     do isol=1,nsol
!!$        psigrid_naked_2(:,isol)=0.0_gp
!!$        print *, "solving for isol ", isol , " ngrid_box_2 ", ngrid_box_2
!!$        call schro(Egrid(isol),rgrid,potgrid,dumgrid1,psigrid_naked_2(1,isol),ngrid_box_2,isol+labs,labs*1.0_gp,zion)
!!$     enddo
!!$     H_2(:,:)=0.0D0
!!$     do i=1,Nsol
!!$        H_2(i,i)=Egrid(i)
!!$        do iproj=1,3
!!$           do igrid=1,Ngrid
!!$              dumgrid1(igrid)=psigrid_naked_2(igrid,i)*projgrid(igrid,iproj)
!!$           enddo
!!$           call integrate(dumgrid1,dumgrid2,rgrid,ngrid_box_2)
!!$           ppgrid(i,iproj)=dumgrid2(ngrid_box_2)
!!$        enddo
!!$     enddo
!!$
!!$     do i=1,Nsol
!!$        do j=1, Nsol
!!$           if ( labs.le.lpx) then
!!$              H_2(i,j)=H_2(i,j)+ ppgrid(i,1)*hsep(1,labs+1)*ppgrid(j,1)&
!!$                   + ppgrid(i,1)*hsep(2,labs+1)*ppgrid(j,2)&
!!$                   + ppgrid(i,2)*hsep(2,labs+1)*ppgrid(j,1)&
!!$                   + ppgrid(i,2)*hsep(3,labs+1)*ppgrid(j,2)&
!!$                   + ppgrid(i,1)*hsep(4,labs+1)*ppgrid(j,3)&
!!$                   + ppgrid(i,3)*hsep(4,labs+1)*ppgrid(j,1)&
!!$                   + ppgrid(i,2)*hsep(5,labs+1)*ppgrid(j,3)&
!!$                   + ppgrid(i,3)*hsep(5,labs+1)*ppgrid(j,2)&
!!$                   + ppgrid(i,3)*hsep(6,labs+1)*ppgrid(j,3)
!!$           endif
!!$           do igrid=1,ngrid_box_2
!!$              dumgrid1(igrid)=psigrid_naked_2(igrid,i)*psigrid_naked_2(igrid,j)*vxcgrid(igrid)
!!$           enddo
!!$           call integrate(dumgrid1,dumgrid2,rgrid,ngrid_box_2)
!!$           H_2(i,j)=H_2(i,j)+dumgrid2(ngrid_box_2)
!!$        enddo
!!$     enddo
!!$     call DSYEV('V','U', Nsol, H_2, Nsol,Egrid_tmp_2 , WORK, Nsol*Nsol*2, INFO)
!!$     call  DGEMM('N','N',Ngrid ,Nsol,Nsol,1.0d0,psigrid_naked_2,Ngrid,&
!!$          H_2 ,Nsol, 0.0D0 , psigrid_not_fitted_2 , Ngrid)
!!$  end if

  if(dofit) then

     print *, "doing the fit "
     Hcorrected=H-Hadd*0.000_gp
     call DSYEV('V','U', Nsol, Hcorrected, Nsol,Egrid_tmp , WORK, Nsol*Nsol*2, INFO)

     real_start=-1
     do iocc=1, Nsol
        !! the conditio below relies on the fact that hgh fitted energies
        !! are good within 10**-3
        if((Etofit(iocc)+0.1).ge.Egrid_tmp(1)) then
           real_start = iocc
           exit
        endif
     enddo

     print *, "real start ", real_start
     Nsol_used=Nsol-( real_start-1  )

     if(.true. )  then
        print *, " routine gatom_modified  ,  comparaison between  energies real and  pseudo-not_fitted "
        do iocc=1, Nsol
           if(iocc.lt.real_start) then
              print *,  iocc, Etofit(iocc) 
           else
              print *,  iocc, Etofit(iocc) , Egrid_tmp(iocc-real_start +1)
           endif
        enddo
     endif
     do isol=1,Nsol-real_start+1
        fact_add=0.0_gp
        do  volta=1,4
           
           Hcorrected=H+Hadd*(-0.001_gp+fact_add)
           call DSYEV('V','U', Nsol, Hcorrected, Nsol,Egrid_tmp , WORK, Nsol*Nsol*2, INFO)
           ene_m= Egrid_tmp(isol)

           Hcorrected=H+Hadd*(0.001_gp+fact_add)
           call DSYEV('V','U', Nsol, Hcorrected, Nsol,Egrid_tmp , WORK, Nsol*Nsol*2, INFO)
           ene_p= Egrid_tmp(isol)
           
           fact_add = fact_add+(Etofit(isol+real_start-1 ) -(ene_p+ene_m)/2.0_gp)/((ene_p-ene_m)/0.002_gp)
           

        enddo
           
        Hcorrected=H+Hadd*fact_add
        call DSYEV('V','U', Nsol, Hcorrected, Nsol,Egrid_tmp , WORK, Nsol*Nsol*2, INFO)
        
        call  DGEMM('N','N',Ngrid,1,Nsol,1.0d0,psigrid_naked,Ngrid,&
             Hcorrected(1,isol),Nsol,0.0D0,psigrid(1,isol),Ngrid)
        Egrid(isol)=Egrid_tmp(isol)
        print *, " Egrid , fact_add " ,  Egrid(isol) , fact_add
        
     enddo


     !! scale psigrid so that it matches AE wavefunctions close to the box border
     do isol=1,Nsol-real_start+1
        fixfact = valuesatp(isol+real_start-1)/ psigrid(Ngrid_box-10,isol)
        do igrid=1, Ngrid
           psigrid(igrid,isol)=psigrid(igrid,isol)*fixfact
        enddo
     enddo



     if(present(psipsigrid)) then
        print *, "copy to psipsigrid  "
        psipsigrid=psigrid
     endif
  

     !! this overlap matrix will be used to get the dual functions (ptildes)
     Soverlap=0.0_gp
     do isol=1,Nsol-real_start+1
        do jsol=isol,Nsol-real_start+1
           do igrid=1,Ngrid_box
              dumgrid1(igrid)=psigrid(igrid,isol)*psigrid(igrid,jsol)
           enddo
           call integrate(dumgrid1,dumgrid2,rgrid,Ngrid_box)
           Soverlap(isol,jsol) = dumgrid2(Ngrid_box)
           Soverlap(jsol,isol) = dumgrid2(Ngrid_box)
        end do
     end do
     do isol=Nsol-real_start+2, Nsol
        !! completes the missing value. Anyway Npaw will be  used to limit the dimension
        Soverlap(isol,isol)=1.0_gp
     enddo
     
     !! now get the dual
     get_duality: if(.true.) then
        !! Hcorrect/Hadd is used here as dummy work array
        !! dumgrid1 dumgrid2 and dumgrid3 too
        !! the latter is already dimensioned
        !! to Ngrid and this dimension must be larger than 3*Nsol
        !! as required by DPOSVX. This should be comfortable
        !! In any case we check
        if(Ngrid<3*Nsol) then
           stop  " Ngrid<3*Nsol for DPOSVX dummies in routine gatom_modified "
        endif
        Tpsigrid=transpose(psigrid)
        !! Nota bene :  we solve only Npaw equations and set preventively the result to zero 
        !! for the other lines of Tpsigrid_dum
        Tpsigrid_dum=0.0_gp
        call DPOSVX( 'N', 'U' , Npaw , Ngrid , Soverlap, Nsol , Hadd, Nsol, EQUED, &
             dumgrid1 , Tpsigrid, Nsol , Tpsigrid_dum, Nsol, RCOND, dumgrid2 , dumgrid3, Hcorrected,&
             IWORK, INFO )
        write(*,'(A,1x,E10.4)') "DUALITY : CONDITION NUMBER FROM DPOSVX ", RCOND
        
        if (INFO.ne.0) then
           stop  " INFO.ne.0 from DPOSVX in routine gatom_modified "
        endif
        
        if(.true.) then
           !! here we check  the duality 
           do isol=1,Npaw
              do jsol=isol,Npaw
                 dumgrid1(:)=Tpsigrid_dum(isol,:)
                 dumgrid2=dumgrid1*psigrid(:,jsol)
                 call integrate(dumgrid2,dumgrid1, rgrid, Ngrid_box)
                 print *, " duality ",isol,jsol, dumgrid1(Ngrid_box) 
              end do
           enddo
        endif
        
        !! The operation below has been postponed after the block below
        !! in order to use temporarily psigrid as auxiliary array
        !! psigrid will be the exit wavefunctions from the routine
        !!psigrid = transpose(Tpsigrid_dum)
        !! ------------------------------
        
        
        !! get the patch
        !! Resolve Htilde to get non-fitted eigenvectors.
        !! In this basis the action of Hnonpatched is given
        !! by the non-fitted eigenenergies.
        !! Get the Matrix M(i,j)=scalar(psinonfit(i), ptilde(j)).
        !! In terms of psinonfit the actions of Htilde on a ptilde(j)
        !! has coefficient MM(i,j)
        !!                Enonfit(i)  M(i,j)
        !!  
        !! The action MMM of Hnonfit in the ptilde basis satisfies
        !!
        !!         M*MMM = MM
        !! This will be solved here with  DPOSVX initialising 
        !! minus PAWpatch in place of MMM. We will then add
        !! the AE energies to the diagonal so that the obtained patch
        !! corresponds effectively to the difference between Hpaw
        !! and Htilde
        
        !!--  get non-fitted eigenvectors
        !!   H is still the non corrected matrix
        Hcorrected=H
        call DSYEV('V','U', Nsol, Hcorrected, Nsol,Egrid_tmp , WORK, Nsol*Nsol*2, INFO)
        call  DGEMM('N','N',Ngrid ,Nsol,Nsol,1.0d0,psigrid_naked,Ngrid,&
             Hcorrected ,Nsol, 0.0D0 , psigrid_not_fitted , Ngrid)
        


        if(.true.) then
           open(unit=22,file='numerov_pseudo_nonfitted.dat')
           do igrid=1, Ngrid
              write(22,'(200(f20.10,1x))') rgrid(igrid), (psigrid_not_fitted(igrid,j ), j=1,Nsol)   
           enddo
           close(unit=22)
        endif
        



        !! calculate overlap between psitilde isol and psigrid_notfitted  jsol 
        !! psigrid  is psitilde ( while ptilde is still stored in Tpsigrid_dum)
        do isol=1,Npaw
           do jsol=1,Nsol_used
              dumgrid2=psigrid(:,isol)*psigrid_not_fitted(:, jsol)
              call integrate(dumgrid2,dumgrid1, rgrid, Ngrid_box)
              Soverlap(isol,jsol) = dumgrid1(Ngrid_box)
           end do
        enddo

        do isol=1,Npaw
           do jsol=1,Nsol_used
              Hcorrected(jsol,isol)= Egrid_tmp(jsol)*Soverlap(isol,jsol)
           end do
        enddo

        call  DGEMM('N','N',Npaw ,Npaw,Nsol_used,1.0d0,Soverlap ,Nsol,&
             Hcorrected ,Nsol, 0.0D0 , PAWpatch , Npaw)
        
        if(INFO .ne. 0) then
           print *, "INFO ", info
           stop "INFO .ne. 0 in  DGESV"
        endif
        do isol=1,Npaw
           do jsol=1,Npaw
              PAWpatch(isol,jsol)=-PAWpatch(isol,jsol)
           enddo
           PAWpatch(isol,isol)=PAWpatch(isol,isol)+Egrid(isol)  
        enddo
        

        if( .true. ) then
           !!! MEGA-check 
           !!    we have still in H the hgh hamiltonian
           !!    written in the psigrid_non_fitted basis
           !!  We are going to patch it with the patch
           !!  We are getting H' = H + Soverlap^T . Pawpatch . Soverlap
           !!   and we get also   S = Identity +   Soverlap^T . (Identity - Spsitildes  ). Soverlap
           !! Where Soverlap(isol,jsol) is overlap between >>Ptilde<< isol and psigrid_notfitted  jsol
           !!         Spsitildes is the overlap matrix between psitilde and psitilde
           !! Then we get the approximated eigenvalues resolving a generalised eigenproblem
           
           
           !! calculate overlap between ptilde isol and psigrid_notfitted  jsol.
           !! Ptilde is still stored in Tpsigrid_dum
           do isol=1,Npaw
              do jsol=1,Nsol_used
                 dumgrid2=Tpsigrid_dum(isol, :)*psigrid_not_fitted(:, jsol)
                 call integrate(dumgrid2,dumgrid1, rgrid, Ngrid_box)
                 Soverlap(isol,jsol) = dumgrid1(Ngrid_box)
              end do
           enddo
           
           call  DGEMM('N','N',Npaw,Nsol_used,Npaw,1.0d0,PAWpatch ,Npaw,&
                Soverlap ,Nsol, 0.0D0 , dumH , Nsol)

           call  DGEMM('T','N',Nsol_used,Nsol_used,Npaw,1.0d0,Soverlap ,Nsol,&
                dumH ,Nsol, 0.0D0 , genH , Nsol)

           do jsol=1,Nsol_used
              genH(jsol,jsol) = genH(jsol,jsol)  +Egrid_tmp(jsol)   !! + H
           end do


           Spsitildes=0.0_gp
           do isol=1,Npaw
              do jsol=isol,Npaw
                 do igrid=1,Ngrid_box
                    dumgrid1(igrid)=psigrid(igrid,isol)*psigrid(igrid,jsol)
                 enddo
                 call integrate(dumgrid1,dumgrid2,rgrid,Ngrid_box)
                 Spsitildes(isol,jsol) = dumgrid2(Ngrid_box)
                 Spsitildes(jsol,isol) = dumgrid2(Ngrid_box)
              end do
              Spsitildes(isol,isol) =Spsitildes(isol,isol) -1.0
           end do



           
           call  DGEMM('N','N',Npaw,Nsol_used,Npaw,1.0d0,Spsitildes ,Npaw,&
                Soverlap ,Nsol, 0.0D0 , dumH , Nsol)
           call  DGEMM('T','N',Nsol_used,Nsol_used,Npaw,1.0d0,Soverlap ,Nsol,&
                dumH ,Nsol, 0.0D0 , genS , Nsol)

           do isol=1,Nsol_used
              do jsol=1,Nsol_used
                 genS(isol,jsol)=-genS(isol,jsol)
              end do
              genS(isol,isol) =1+ genS(isol,isol)
           end do
           
    
!!$           do isol=1,Npaw
!!$              do jsol=1,Nsol_used
!!$                 dumgrid2=psigrid(:, isol)*psigrid_not_fitted(:, jsol)
!!$                 call integrate(dumgrid2,dumgrid1, rgrid, Ngrid_box)
!!$                 Soverlap(isol,jsol) = dumgrid1(Ngrid_box)
!!$              end do
!!$           enddo
!!$
!!$
!!$           call  DGEMM('N','T',Nsol_used,Npaw,Nsol_used,1.0d0,genH ,Nsol,&
!!$                Soverlap ,Nsol, 0.0D0 , dumH , Nsol)
!!$
!!$           call  DGEMM('N','N',Npaw,Npaw,Nsol_used,1.0d0,Soverlap ,Nsol,&
!!$                dumH ,Nsol, 0.0D0 , genH , Nsol)
!!$
!!$           do i=1,Npaw
!!$              print *, genH(i,i)
!!$           end do
!!$           do i=1,Npaw
!!$              print *, (genH(i,j),j=1,Npaw)
!!$           end do
!!$           
!!$           stop

!!$
!!$
!!$           call  DGEMM('N','T',Nsol_used,Npaw,Nsol_used,1.0d0,genS ,Nsol,&
!!$                Soverlap ,Nsol, 0.0D0 , dumH , Nsol)
!!$
!!$           call  DGEMM('N','N',Npaw,Npaw,Nsol_used,1.0d0,Soverlap ,Nsol,&
!!$                dumH ,Nsol, 0.0D0 , genS , Nsol)
!!$
!!$           do i=1,Npaw
!!$              print *, genS(i,i)
!!$           end do
!!$           do i=1,Npaw
!!$              print *, (genS(i,j),j=1,Npaw)
!!$           end do
!!$           
!!$
!!$
!!$           stop

           ITYPE=1
           LDWORK=Ngrid
           CALL  DSYGV(ITYPE, "N", "U", Nsol_used, genH, Nsol, genS, Nsol , dumgrid2 , dumgrid3, &
                LDWORK, INFO)
           print *, " first  eigenvalues "
           do i=1,Npaw
              print *, dumgrid2( i)
           end do
        endif


!!$        if( present(psipsigrid) ) then
!!$
!!$           do isol=1,Npaw
!!$              do jsol=1,Nsol_used
!!$                 dumgrid2=Tpsigrid_dum(isol, :)*psigrid_not_fitted_2(:, jsol)
!!$                 call integrate(dumgrid2,dumgrid1, rgrid, Ngrid_box)
!!$                 Soverlap(isol,jsol) = dumgrid1(Ngrid_box)
!!$              end do
!!$           enddo
!!$           
!!$           call  DGEMM('N','N',Npaw,Nsol_used,Npaw,1.0d0,PAWpatch ,Npaw,&
!!$                Soverlap ,Nsol, 0.0D0 , dumH , Nsol)
!!$
!!$           call  DGEMM('T','N',Nsol_used,Nsol_used,Npaw,1.0d0,Soverlap ,Nsol,&
!!$                dumH ,Nsol, 0.0D0 , genH , Nsol)
!!$
!!$           do jsol=1,Nsol_used
!!$              genH(jsol,jsol) = genH(jsol,jsol)  +Egrid_tmp_2(jsol)   !! + H
!!$           end do
!!$           
!!$           call  DGEMM('N','N',Npaw,Nsol_used,Npaw,1.0d0,Spsitildes ,Npaw,&
!!$                Soverlap ,Nsol, 0.0D0 , dumH , Nsol)
!!$           call  DGEMM('T','N',Nsol_used,Nsol_used,Npaw,1.0d0,Soverlap ,Nsol,&
!!$                dumH ,Nsol, 0.0D0 , genS , Nsol)
!!$           do isol=1,Nsol_used
!!$              do jsol=1,Nsol_used
!!$                 genS(isol,jsol)=-genS(isol,jsol)
!!$              end do
!!$              genS(isol,isol) =1+ genS(isol,isol)
!!$           end do
!!$
!!$           ITYPE=1
!!$           LDWORK=Ngrid
!!$           CALL  DSYGV(ITYPE, "V", "U", Nsol_used, genH, Nsol, genS, Nsol , dumgrid2 , dumgrid3, &
!!$                LDWORK, INFO)
!!$           print *, " first  eigenvalues , INFO", INFO
!!$
!!$           do i=1,Npaw
!!$              print *, dumgrid2( i)
!!$           end do
!!$
!!$           call  DGEMM('N','N',Ngrid, Nsol_used,  Nsol_used  ,1.0d0,  psigrid_not_fitted_2    ,Ngrid,&
!!$              genH   ,Nsol, 0.0D0 , psipsigrid ,Ngrid )
!!$        endif


        !! these are the wavefunctions ptilde returned by the routine
        psigrid = transpose(Tpsigrid_dum)
        
     endif get_duality
  else
     call DSYEV('V','U', Nsol, H, Nsol,Egrid , WORK, Nsol*Nsol*2, INFO)
     call  DGEMM('N','N',Ngrid ,Nsol,   Nsol,1.0d0 ,psigrid_naked, Ngrid ,H,Nsol, 0.0D0 , psigrid , Ngrid)
  endif



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
!!!        write(66,*)  lmax+1
!!!        write(66,*) ' #LINETYPE{1324}' 
!!!        write(66,*) ' $' 
!!!  do l=0,lmax
!!!           write(66,*) ' 161'
!!!     r=0._gp
!!!     do
!!!        tt= wave(ng,l,xp,psi(0,1,l+1),r)
!!!              write(66,*) r,tt
!!!        r=r+.025_gp
!!!        if(r > 4.00001_gp) exit
!!!     end do
!!!  end do
! writing lines suppressed
!!!        write(67,*) min(lmax+1,3)
!!!        write(67,*) ' #LINETYPE{132}'
!!!        write(67,*) ' #TITLE{FOURIER}' 
!!!        write(67,*) ' $'
  dr=6.28_gp/rprb/200._gp
!!!        write(67,*) ' 200'
  rk=0._gp
  loop_rk1: do 
     tt=0._gp
     do i=0,ng
        texp=exp(-.25_gp*rk**2/xp(i))
!        texp=exp(-.5_gp*energy/xp(i))
        sd=sqrt(xp(i))
        tt=tt+psi(i,1,1)*0.4431134627263791_gp*texp/sd**3
     end do
!!!           write(67,*) rk,tt
     rk=rk+dr
     if(rk > 6.28_gp/rprb-.5_gp*dr) exit loop_rk1
  end do loop_rk1
  if (lmax.ge.1) then
!!!           write(67,*) ' 200'
     rk=0._gp
     loop_rk2: do 
        tt=0._gp
        do i=0,ng
           texp=exp(-.25_gp*rk**2/xp(i))
           sd=sqrt(xp(i))
           tt=tt+psi(i,1,2)*0.2215567313631895_gp*rk*texp/sd**5
        end do
!!!              write(67,*) rk,tt
        rk=rk+dr
        if (rk > 6.28_gp/rprb-.5_gp*dr) exit loop_rk2
     end do loop_rk2
  end if
  if (lmax.ge.2) then
!!!           write(67,*) ' 200'
     rk=0._gp
     do 
        tt=0._gp
        do i=0,ng
           texp=exp(-.25_gp*rk**2/xp(i))
           sd=sqrt(xp(i))
           tt=tt+psi(i,1,3)*0.1107783656815948_gp*rk**2*texp/sd**7
        end do
!!!              write(67,*) rk,tt
        rk=rk+dr
        if (rk > 6.28_gp/rprb-.5_gp*dr) exit
     end do
  end if
END SUBROUTINE gatom_modified


!! this routine solves exactly a reference model
subroutine gatom_modified_eqdiff(rcov,rprb,lmax,lpx,noccmax,occup,&
                 zion,alpz,gpot,alpl,hsep,alps,vh,xp,rmt,fact,nintp,&
                 aeval,ng,psi,res,chrg,&
                 Nsol, Labs, Ngrid,Egrid,  rgrid , psigrid )
  use module_base
  use esatto

  implicit real(gp) (a-h,o-z)
  logical :: noproj, readytoexit
  integer, parameter :: n_int=1000
  dimension psi(0:ng,noccmax,lmax+1),aeval(noccmax,lmax+1),&
       hh(0:ng,0:ng),ss(0:ng,0:ng),eval(0:ng),evec(0:ng,0:ng),&
       gpot(3),hsep(6,lpx+1),rmt(n_int,0:ng,0:ng,lmax+1),&
       pp1(0:ng,lpx+1),pp2(0:ng,lpx+1),pp3(0:ng,lpx+1),alps(lpx+1),&
       potgrd(n_int),&
       rho(0:ng,0:ng,lmax+1),rhoold(0:ng,0:ng,lmax+1),xcgrd(n_int),&
       occup(noccmax,lmax+1),chrg(noccmax,lmax+1),&
       vh(0:ng,0:ng,4,0:ng,0:ng,4),&
       res(noccmax,lmax+1),xp(0:ng),& 
       psigrid(Ngrid, Nsol),psigrid_naked(Ngrid,Nsol), projgrid(Ngrid,3), &
       rhogrid(Ngrid), potgrid(Ngrid),  newpotgrid(Ngrid), &
       vxcgrid(Ngrid), &
       Egrid(nsol), ppgrid(Nsol,3), work(nsol*nsol*2), &
       H(Nsol, Nsol)
  
  real(gp) Ediff, Rdiff, rpot_a,spot_a,hpot_a, Rinf_a

  real(gp) :: rgrid(Ngrid)
  real(gp), target :: dumgrid1(Ngrid),dumgrid2(Ngrid) ,dumgrid3(Ngrid)
  real(gp) , pointer :: y_r(:) ,  d_r(:)
  character(len=*), parameter :: subname='gatom_modified_eqdiff'

 

!!!
!!!  y_r(0)=1.0
!!!  y_r(1)=2
!!!  y_r(2)=3
!!!  !! py_r=> y_r
!!!  print *, py_r(1)
!!! 

  if (nintp.ne.n_int) stop 'n_int><nintp  xabs'


  do l=0,lmax
     if (occup(1,l+1).gt.0._gp) lcx=l
  end do
  !write(6,*) 'lcx',lcx
 
  noproj=.true.
  do l=1,lpx+1
     noproj = noproj .and. (alps(l) .eq. 0._gp)
  end do
  
  noproj=.true.

! projectors, just in case
  if ( .not. noproj) then
     do l=0,lpx
        gml1=sqrt( gamma_restricted(real(l,gp)+1.5_gp) / (2._gp*alps(l+1)**(2*l+3)) )
        gml2=sqrt( gamma_restricted(real(l,gp)+3.5_gp) / (2._gp*alps(l+1)**(2*l+7)) )&
            /(real(l,gp)+2.5_gp)
        gml3=sqrt( gamma_restricted(real(l,gp)+5.5_gp) / (2._gp*alps(l+1)**(2*l+11)) )&
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
                 projgrid(igrid,iorder) = (  exp( -r*r*tt)*(r**l *r) )* r**(2*(iorder-1))
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
  readytoexit=.false.
  big_loop: do it=1,50
     if(it.eq.50) readytoexit=.true.
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


     if( readytoexit) then

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
     endif
  
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


     if(readytoexit) then
        
        do igrid=1, Ngrid
           r=rgrid(igrid)
           potgrid(igrid) =0.5_gp*r*r  /    rprb**4*0 
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
     endif
! ------------------------------------------------------------------------------------



 

     loop_l: do l=0,lmax
        gml=.5_gp*gamma_restricted(.5_gp+real(l,gp))

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
                   .5_gp*const*sxp**2*(real(l,gp)+.5_gp)*(real(l,gp)+1.5_gp)/rprb**4*0 
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
              hh(i,j)=hh(i,j)+ gpot(1)*.5_gp*gamma_restricted(1.5_gp+real(l,gp))*tt**(1.5_gp+real(l,gp))&
                   + (gpot(2)/alpl**2)*.5_gp*gamma_restricted(2.5_gp+real(l,gp))*tt**(2.5_gp+real(l,gp))&
                   + (gpot(3)/alpl**4)*.5_gp*gamma_restricted(3.5_gp+real(l,gp))*tt**(3.5_gp+real(l,gp))
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
     ! write(6,*) 'evdiff',it,tt, readytoexit
     if (tt.lt.1.e-12_gp) then
        if( readytoexit) then
           exit big_loop
        endif
        readytoexit=.true.
     end if
  end do big_loop
! End of the big loop




  

  if(.true.) then
     nls_a=14
     allocate(y_r(0:nls_a-1+ndebug), stat=i_stat)
     call memocc(i_stat,y_r,'y_r',subname)
     allocate(d_r(0:nls_a-1+ndebug), stat=i_stat)
     call memocc(i_stat,d_r,'d_r',subname)
 

     if(.false.) then ! ricalcola la rho con le soluzioni radiali

        do iter=1,50
           
           do igrid=1,Ngrid
              rhogrid(igrid)=0.0_gp
           enddo
           do l=0, lmax
              dumgrid3 = potgrid  +   0.5_gp*l *(l +1.0_gp)/rgrid/rgrid  + vxcgrid
              
              do igrid=1,Ngrid

                 if(-1.0/rgrid(igrid)< potgrid(igrid)  ) dumgrid3(igrid)=&
                      -1.0/rgrid(igrid)   +  &
                      0.5_gp*l *(l +1.0_gp)/rgrid(igrid)/rgrid(igrid)  + vxcgrid(igrid)
                 !!! if(rgrid(igrid)>8.0+iter/20) dumgrid3(igrid)=10.0
              enddo
              
              dumgrid1(:)=0.0_gp
              do iocc=1,noccmax
                 if( occup(iocc,l+1)>0.0001 )  then
                    call schro(Egrid(isol) , rgrid ,  dumgrid3 , dumgrid1, dumgrid2 , ngrid , iocc+l , l*1.0_gp ,  zion)

                    do igrid=1,Ngrid
                       dum=dumgrid2(igrid)
                       rhogrid(igrid)=rhogrid(igrid)+dum*dum*occup(iocc,l+1)
                    enddo
                 endif
              enddo
           enddo


           do igrid=1,Ngrid
              r=rgrid(igrid)
              dum = rhogrid(igrid)/r/r *0.07957747154594768_gp
              vxcgrid(igrid)=0.8_gp*vxcgrid(igrid) +0.2_gp* emuxc(dum) 
           enddo

           do igrid=1, Ngrid
              r=rgrid(igrid)
              newpotgrid(igrid) =0.5_gp*r*r  /    rprb**4*0 
              newpotgrid(igrid) = newpotgrid(igrid) - zion/r * derf( r/alpz/sqrt(2.0)   )
              rr = r/alpz
              newpotgrid(igrid) = newpotgrid(igrid) + exp(-0.5_gp * rr**2 )*( gpot(1)+gpot(2)*rr**2 + gpot(3)*rr**4 )
           enddo
           ! poisson per potgrid
           call integrate( rhogrid, dumgrid1, rgrid, Ngrid)
           do igrid=1, Ngrid
              newpotgrid(igrid)=newpotgrid(igrid)+dumgrid1(igrid)/rgrid(igrid)
           enddo

           do igrid=1, Ngrid
              dumgrid1(igrid) = rhogrid(igrid)/rgrid(igrid)
           enddo
           call integrate( dumgrid1, dumgrid2, rgrid, Ngrid)
           do igrid=1, Ngrid
              newpotgrid(igrid)=newpotgrid(igrid)-dumgrid2(igrid) + dumgrid2(Ngrid)
           enddo

           potgrid=0.8_gp*potgrid+0.2_gp*newpotgrid

        enddo
     endif
     
     dumgrid1(:)=0.0_gp
     open(unit=22,file='pot2.dat')
     do igrid=1, Ngrid
        write(22,'(200(f20.10,1x))') rgrid(igrid),potgrid(igrid)+vxcgrid(igrid)
     enddo
     close(unit=22)
     
     print *, " trying  eq diff "
     do idiff=1,600
        Ediff=idiff*0.04_gp/3.0 +0.0_gp
        do igrid=Ngrid,1,-1
           dumgrid2(igrid)=potgrid(igrid)+vxcgrid(igrid)
           if(rgrid(igrid)>2.0_gp) then
              Nrdiff = igrid;
              Rdiff = rgrid(igrid)
           endif
        enddo
        do l=0, nls_a-1
           if(l>0) then
              do igrid=1,Ngrid
                 dumgrid2 (igrid)=dumgrid2 (igrid)+ 0.5_gp*(2.0_gp*l)/rgrid(igrid)/rgrid(igrid)
              enddo
           endif
           d_r(l)=phase( Ediff, Nrdiff, rgrid,dumgrid2  , dumgrid1 , psigrid(1,2+l) , l*1.0_gp ,0, 1)
           y_r(l)= psigrid(Nrdiff,2+l)
           
        enddo
        
        dumgrid1 = 0.0_gp
        dumgrid2 = potgrid + vxcgrid
        call schro(Egrid(1) , rgrid , dumgrid2  , dumgrid1, psigrid(:,1) , ngrid , 1+0 , 0*1.0_gp ,  zion)
        
        open(unit=22,file='bound2.dat')
        do igrid=1, Ngrid
           write(22,'(200(f20.10,1x))') rgrid(igrid),  psigrid(igrid,1), psigrid(igrid,2+Labs)
        enddo
        close(unit=22)
        
        dumgrid1(:) =  psigrid(:,2+Labs)* psigrid(:,1)
        
        do l=1,Labs
           dumgrid1(:) = dumgrid1(:) *rgrid(:)
        end do
        call integrate(dumgrid1(1),dumgrid3(1),rgrid(1) ,Nrdiff)
        fattore = dumgrid3(Nrdiff)

        dumgrid1(:)=0
        if(.true.) then
           lpot_a=1
           lpot_a=1
           rpot_a = 6.0_gp
           spot_a = 1.0_gp
           hpot_a = 3.0_gp
           Rinf_a=100.0_gp
           nstesp_coarse=1000
           nsteps_fine  = 40
                     
!!!           print *, y_r,d_r
!!!           print *, " A  y_r(1) " , y_r(1), d_r(1)

           ref= esatto_CalcolaRiflettivita( ngrid, rgrid, dumgrid2, nls_a, lpot_a, rpot_a,spot_a,hpot_a,y_r,d_r,&
                Rdiff,    Rinf_a ,nstesp_coarse ,nsteps_fine, Ediff , Labs)

       

           print *," ",  Ediff, ref*fattore*fattore, fattore
           

!!!           do igrid=1,Ngrid
!!!              dumgrid2 (igrid)=potgrid(igrid)+vxcgrid(igrid)+ 0.5_gp*(2.0_gp)/rgrid(igrid)/rgrid(igrid)
!!!           enddo
!!!           d_r(1)=phase( Ediff, Ngrid, rgrid,dumgrid2  , dumgrid1 , psigrid(1,2+1) , 1*1.0_gp ,0, 1)
!!!           y_r(1)= psigrid(Ngrid,2+1)
!!!
!!!           print *,"C2 ", Ediff, 4*fattore*fattore/(   y_r(1)*y_r(1) +  d_r(1)*d_r(1)/2.0/Ediff   )


!!!           open(unit=22,file='libera.dat')
!!!           do igrid=1, Ngrid
!!!              write(22,'(200(f20.10,1x))') rgrid(igrid),  psigrid(igrid,2+1)
!!!           enddo
!!!           close(unit=22)
!!!           stop


        endif
     enddo
     i_all=-product(shape(y_r))*kind(y_r)
     deallocate(y_r, stat=i_stat)
     call memocc(i_stat,i_all,'y_r',subname)

     i_all=-product(shape(d_r))*kind(d_r)
     deallocate(d_r, stat=i_stat)
     call memocc(i_stat,i_all,'d_r',subname)
     stop
  endif



  do igrid=1, ngrid
     r=rgrid(igrid)
     potgrid(igrid)=potgrid(igrid)+ 0.5_gp*labs*(labs+1.0_gp)/r/r
  enddo
  
  dumgrid1(:)=0.0_gp
  do isol=1,nsol
     call schro(Egrid(isol) , rgrid ,  potgrid , dumgrid1, psigrid_naked(1,isol) , ngrid, isol+labs , labs*1.0_gp ,  zion)

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
!!!        write(66,*)  lmax+1
!!!        write(66,*) ' #LINETYPE{1324}' 
!!!        write(66,*) ' $' 
!!!  do l=0,lmax
!!!           write(66,*) ' 161'
!!!     r=0._gp
!!!     do
!!!        tt= wave(ng,l,xp,psi(0,1,l+1),r)
!!!              write(66,*) r,tt
!!!        r=r+.025_gp
!!!        if(r > 4.00001_gp) exit
!!!     end do
!!!  end do
! writing lines suppressed
!!!        write(67,*) min(lmax+1,3)
!!!        write(67,*) ' #LINETYPE{132}'
!!!        write(67,*) ' #TITLE{FOURIER}' 
!!!        write(67,*) ' $'
  dr=6.28_gp/rprb/200._gp
!!!        write(67,*) ' 200'
  rk=0._gp
  loop_rk1: do 
     tt=0._gp
     do i=0,ng
        texp=exp(-.25_gp*rk**2/xp(i))
!        texp=exp(-.5_gp*energy/xp(i))
        sd=sqrt(xp(i))
        tt=tt+psi(i,1,1)*0.4431134627263791_gp*texp/sd**3
     end do
!!!           write(67,*) rk,tt
     rk=rk+dr
     if(rk > 6.28_gp/rprb-.5_gp*dr) exit loop_rk1
  end do loop_rk1
  if (lmax.ge.1) then
!!!           write(67,*) ' 200'
     rk=0._gp
     loop_rk2: do 
        tt=0._gp
        do i=0,ng
           texp=exp(-.25_gp*rk**2/xp(i))
           sd=sqrt(xp(i))
           tt=tt+psi(i,1,2)*0.2215567313631895_gp*rk*texp/sd**5
        end do
!!!              write(67,*) rk,tt
        rk=rk+dr
        if (rk > 6.28_gp/rprb-.5_gp*dr) exit loop_rk2
     end do loop_rk2
  end if
  if (lmax.ge.2) then
!!!           write(67,*) ' 200'
     rk=0._gp
     do 
        tt=0._gp
        do i=0,ng
           texp=exp(-.25_gp*rk**2/xp(i))
           sd=sqrt(xp(i))
           tt=tt+psi(i,1,3)*0.1107783656815948_gp*rk**2*texp/sd**7
        end do
!!!              write(67,*) rk,tt
        rk=rk+dr
        if (rk > 6.28_gp/rprb-.5_gp*dr) exit
     end do
  end if



END SUBROUTINE gatom_modified_eqdiff






subroutine GetExcitedOrbitalAsG( in_iat_absorber , atoms, rxyz, nproc, iproc,&
     abs_initial_N, abs_initial_L, abs_final_L,abs_rpower,in_NPaw, do_eqdiff, &
     filename, ha )

  use module_base
  use module_types
  use module_interfaces

  implicit none
  integer :: abs_final_L 
  integer, intent(in) :: in_iat_absorber, nproc, iproc
  !!$ type(gaussian_basis) , intent(out) :: Gabsorber
  type(atoms_data), intent(in) :: atoms
  real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
  !! real(wp), dimension(2*abs_final_L+1), intent(out) :: Gabs_coeffs
  logical do_eqdiff
  character(len=*) filename
  character(len=1000) filename_in_coeffs
  type(lanczos_args) :: ha
  integer ::  abs_initial_L, abs_initial_N, abs_rpower , in_NPaw
  ! -----------------------------------------------------------
  
  real(gp) :: coeffs_in_paw(in_Npaw)
  integer :: ity, ng , noccmax, lmax, j, ierr, i_all
  real(gp) , pointer :: expo(:), psi(:,:,:),psireal(:,:,:), aeval(:,:), occup(:,:), gcoeffs(:)
  integer :: psp_modifier
  integer :: ig, iocc, iexpo


  integer, parameter :: Norder=4, dump_functions=1

  real(gp) , pointer:: psi1s(:) 
  integer :: real_start
  real(gp) :: cradius, boxradius
  
  integer ::  Nsol , Ngrid, igrid, Ngrid_box,   Npaw    ,Ngrid_biggerbox 

  integer :: ng_fine
  real(gp), pointer :: expo_fine(:)

  integer i_stat
  character(len=*), parameter :: subname='GetExcitedOrbitalAsG'
  real(gp) rzero
  real(gp), pointer ::  Egrid(:) , Egrid_pseudo(:) 
  real(gp), pointer::PAWpatch(:,:),rgrid(:),psigrid(:,:),psigrid_pseudo(:,:) &
       ,psipsigrid_pseudo(:,:), psigrid_bigger(:,:)
  integer isol, jsol
  character(len=2) :: symbol
  real(gp), dimension(6,4) :: neleconf
  integer :: nsccode,mxpl,mxchg
  real(kind=8) ::rcov,rprb,ehomo, amu
  real(kind=8) sum



  ! if (in_iat_absorber.ne.0) then

  ity = atoms%iatype(in_iat_absorber)
  !! ATTENTION put this to 60 when doing exact model with gatom ( instead of atom-numeric)
  ng  = 30
  noccmax = 5 
  lmax=3
  
  ng_fine= 350
  
  Nsol=200

  !! if Npaw !=Nsol
  !! the dual functions, PAWpatch will be calculated,
  !!  otherwise old style calculation will be done
  !! Npaw=Nsol

  if(in_NPaw.eq.0) then
     Npaw=Nsol
  else
     Npaw=in_NPaw
  endif
  
  Ngrid=20000
  
  if(do_eqdiff) then
     cradius=30.0_gp !!!!!!!! ATTENZIONE
     Ngrid=60000
  else
     cradius=10.0_gp !!!!!!!! ATTENZIONE
  endif
  !!!!  the cradius is used to calculate self consistent atom numerically
  !!!! and solutions of the numerical equation pseudo and real.
  !!!! When searching for the projector basis, however
  !!!!  the effective radius is effectively reduced to a minor extent
  !!!! by using a part of the grid only : from 1 up to Ngrid_box ( see below)
  !!!! Beyond Ngrid_box all phi_tilde and p_tilde will be zero
  !!!! Ngrid_box is set approximatively to boxradius
  !! boxradius=3.0


  call eleconf(atoms%nzatom(ity) ,atoms%nelpsp(ity) ,symbol,rcov,rprb,ehomo,&
       neleconf,nsccode,mxpl,mxchg,amu)
  boxradius=rcov*1.0_gp

  print *, "boxradius has been set to ", boxradius
  
  allocate(expo_fine(ng_fine  +ndebug ), stat=i_stat)
  call memocc(i_stat,expo_fine,'expo_fine',subname)
  
  allocate(expo(ng +ndebug  ), stat=i_stat)
  call memocc(i_stat,expo,'expo',subname)
  
  allocate(psi ( 0:ng-1  ,noccmax,lmax+1+ndebug ), stat=i_stat)
  call memocc(i_stat,psi,'psi',subname)
  
  allocate(psireal ( 0:Ngrid  ,noccmax,lmax+1+ndebug ), stat=i_stat)
  call memocc(i_stat,psireal,'psireal',subname)
  
  allocate(gcoeffs ( 0:ng_fine-1  +ndebug ), stat=i_stat)
  call memocc(i_stat,gcoeffs,'gcoeffs',subname)
  
  allocate(aeval ( noccmax  ,lmax+1+ndebug ), stat=i_stat)
  call memocc(i_stat,aeval,'aeval',subname)
  
  allocate(occup ( noccmax  ,lmax+1+ndebug ), stat=i_stat)
  call memocc(i_stat,occup,'occup',subname)
  
  allocate( Egrid(Nsol +ndebug ), stat=i_stat)
  call memocc(i_stat,Egrid,'Egrid',subname)
  
  allocate( PAWpatch(Npaw,Npaw +ndebug ), stat=i_stat)
  call memocc(i_stat,PAWpatch,'PAWpatch',subname)
  



  allocate( rgrid(Ngrid +ndebug ), stat=i_stat)
  call memocc(i_stat,rgrid,'rgrid',subname)
 

 
  allocate( psigrid(Ngrid  , Nsol +ndebug ), stat=i_stat)
  call memocc(i_stat,psigrid,'psigrid',subname)
  
  allocate( Egrid_pseudo(Nsol +ndebug ), stat=i_stat)
  call memocc(i_stat,Egrid_pseudo,'Egrid_pseudo',subname)


  
  allocate( psigrid_pseudo(Ngrid  , Nsol +ndebug), stat=i_stat)
  call memocc(i_stat,psigrid_pseudo,'psigrid_pseudo',subname)

  
  allocate( psipsigrid_pseudo(Ngrid  , Nsol +ndebug), stat=i_stat)
  call memocc(i_stat,psipsigrid_pseudo,'psipsigrid_pseudo',subname)

  allocate( psigrid_bigger(Ngrid  , Nsol +ndebug), stat=i_stat)
  call memocc(i_stat,psigrid_bigger,'psigrid_bigger',subname)


  allocate(psi1s( Ngrid +ndebug ), stat=i_stat)
  call memocc(i_stat,psi1s,'psi1s',subname)
   
!!!  do igrid=1, Ngrid
!!!     rgrid(igrid) = igrid*1.0_gp/Ngrid * cradius
!!!  enddo

  rzero = 1.0D-5/Ngrid * cradius 
  do igrid=1, Ngrid
     rgrid(igrid) = rzero*  exp( igrid*   1.0_gp/Ngrid * log( cradius/rzero ))
     if ( rgrid(igrid)<boxradius) then
        Ngrid_box=igrid
     endif
     !! this 1.5 factor is the same as in file init.f90, routine  fillPawProjOnTheFly
     if ( rgrid(igrid)<boxradius*1.5_gp) then
        Ngrid_biggerbox=igrid
     endif
  enddo
  
  
!!! 
!!!  
!!!
!!!
!!!  if(iproc.eq.0)   print * , " routine GetExcitedOrbitalAsG  , calculate pseudo  noccmax, lmax, per verificare  grid  ", noccmax, lmax 
!!!  psp_modifier=0
!!!  call abs_generator_modified(iproc,atoms%nzatom(ity), atoms%nelpsp(ity),atoms%psppar(0,0,ity),&
!!!       atoms%npspcode(ity),ng-1 ,noccmax , lmax , expo,psi,aeval, occup , psp_modifier , &
!!!       Nsol, abs_initial_L , Ngrid,Egrid_pseudo,  rgrid , psigrid_pseudo  )
!!!
!!!
!!!
!!!  if(iproc.eq.0) then
!!!     do iocc=1,2
!!!        print *, " pseudo  Egau, pseudo  Egrid ", aeval(iocc, abs_initial_L+1  ), Egrid_pseudo(iocc)
!!!     enddo
!!!  endif
!!!



  if(do_eqdiff) then
     if(iproc.eq.0)   print * , " routine GetExcitedOrbitalAsG  solving differential equation "
     psp_modifier=-1

     if (.true.) then

        call abs_generator_modified(iproc,atoms%nzatom(ity), atoms%nelpsp(ity),atoms%psppar(:,:,ity),&
             atoms%npspcode(ity),Ngrid ,noccmax , lmax , expo,psireal ,aeval, occup , psp_modifier  , &
             Nsol, abs_final_L , Ngrid,Ngrid_box, Egrid,  rgrid , psigrid , Npaw, PAWpatch )
        
     else

        call abs_generator_modified(iproc,atoms%nzatom(ity), atoms%nelpsp(ity),atoms%psppar(:,:,ity),&
             atoms%npspcode(ity),Ng ,noccmax , lmax , expo,psi,aeval, occup , psp_modifier , &
             Nsol, abs_final_L , Ngrid,Ngrid_box, Egrid_pseudo,  rgrid , psigrid_pseudo,Npaw,  PAWpatch  )


     endif
     stop
  endif
  
!   !! --------------------------------------------------------------------------------------------------------
!   if(iproc.eq.0)   print * , " routine GetExcitedOrbitalAsG  , calculate pseudo  noccmax, lmax, ", noccmax, lmax 
!   psp_modifier=0
!   call abs_generator_modified(iproc,atoms%nzatom(ity), atoms%nelpsp(ity),atoms%psppar(0,0,ity),&
!        atoms%npspcode(ity),ng-1 ,noccmax , lmax , expo,psi,aeval, occup , psp_modifier , &
!        Nsol, abs_final_L , Ngrid, Ngrid_box,Egrid_pseudo,  rgrid , psigrid_pseudo  )
!   if(dump_functions==1)  then
!      if(iproc.eq.0) then
!         do iocc=1,2
!            print *, " pseudo  Egau, pseudo  Egrid ", aeval(iocc, abs_final_L+1  ), Egrid_pseudo(iocc)
!         enddo
!      endif
!   endif
!    if(iproc.eq.0)   print * , " uscito routine GetExcitedOrbitalAsG  , calculate pseudo " 
!   if(iproc.eq.0 .and. dump_functions.eq.1) then
!      if(psp_modifier.eq.0) then
!         call dump_gauwf_on_radgrid("pseudo_wf_radgrid", ng-1,noccmax,lmax,expo,psi)
!      else
!         call dump_gauwf_on_radgrid("real_wf_radgrid",   ng-1,noccmax,lmax,expo,psi)
!      endif
!   endif
!   if(iproc.eq.0 .and. dump_functions.eq.1) then
!      open(unit=22,file='numerov_pseudo.dat')
!      do igrid=1, Ngrid
!         write(22,'(200(f20.10,1x))') rgrid(igrid), (psigrid_pseudo(igrid,j ), j=1,Nsol)   
!      enddo
!      close(unit=22)
!   endif
!   !! ----------------------------------------------------------------------------------------------

  if(iproc.eq.0) print * , " routine GetExcitedOrbitalAsG  , generate  atom to  extract 1S " 
  psp_modifier=1
  !!  atoms%nzatom(ity)+1
  call abs_generator_modified(iproc,atoms%nzatom(ity), atoms%nelpsp(ity),atoms%psppar(:,:,ity),&
       atoms%npspcode(ity),Ngrid ,noccmax , lmax , expo,psireal ,aeval, occup , psp_modifier  , &
       7, abs_initial_L , Ngrid,Ngrid_box, Egrid,  rgrid , psigrid , Npaw, PAWpatch )
  
  !! retrieve 1 s *r 
  do igrid=1,Ngrid
     psi1s(igrid) =   psigrid(igrid,abs_initial_N- abs_initial_L  )  *( rgrid(igrid) ** abs_rpower  )
  enddo

!!$  
!!$  sum=0.0
!!$  do igrid=1,Ngrid-1
!!$     sum=sum+ (rgrid(igrid+1)-rgrid(igrid)) * psigrid(igrid,abs_initial_N- abs_initial_L  )**2
!!$  enddo
!!$  print *, "norma dello stato iniziale ", sum
!!$  sum=0.0
!!$  do igrid=1,Ngrid-1
!!$     sum=sum+ ((rgrid(igrid+1)-rgrid(igrid)) * psigrid(igrid,abs_initial_N- abs_initial_L  )**2)   *&
!!$          (( rgrid(igrid) ** abs_rpower  )**2)
!!$  enddo
!!$  print *, "norma dello stato iniziale con operatore ", sum
 

  if(iproc.eq.0 .and. dump_functions.eq.1) then
     open(unit=22,file='numerov0.dat')
     do igrid=1, Ngrid
        write(22,'(200(f20.10,1x))') rgrid(igrid), (psigrid(igrid,j ), j=1,Nsol)   
     enddo
     close(unit=22)
  endif


  if(iproc.eq.0) print * , " routine GetExcitedOrbitalAsG  , calculate  atom with  real-pot, abs_final_L " ,abs_final_L
  psp_modifier=1
  call abs_generator_modified(iproc,atoms%nzatom(ity), atoms%nelpsp(ity),atoms%psppar(:,:,ity),&
       atoms%npspcode(ity),Ngrid ,noccmax , lmax , expo,psireal,aeval, occup , psp_modifier ,  &
       Nsol, abs_final_L , Ngrid, Ngrid_box,Egrid,  rgrid , psigrid ,Npaw,  PAWpatch )


!!$  if(dump_functions==1)  then
!!$     if(iproc.eq.0) then
!!$        do iocc=1,3
!!$           print *, " real Egau, real Egrid ", aeval(iocc, abs_final_L+1), Egrid(iocc)
!!$        enddo
!!$     endif
!!$  endif


  if(iproc.eq.0 .and. dump_functions.eq.1) then
     open(unit=22,file='numerov.dat')
     do igrid=1, Ngrid
        write(22,'(200(f20.10,1x))') rgrid(igrid), (psigrid(igrid,j ), j=1,Nsol)   
     enddo
     close(unit=22)
  endif

  if(iproc.eq.0 .and. dump_functions.eq.1) then
     if(psp_modifier.eq.0) then
        call dump_gauwf_on_radgrid("pseudo_wf_radgrid", ng-1,noccmax,lmax,expo,psi)
     else
        call dump_real_on_radgrid("real_wf_radgrid",   Ngrid,noccmax,lmax,expo,psireal, rgrid)
     endif
  endif
  
  if(iproc.eq.0)   print * , " routine GetExcitedOrbitalAsG  , calculate pseudo  noccmax, lmax, ", noccmax, lmax 
  psp_modifier=0
  Egrid_pseudo(:)= Egrid(:)  !! to fit these energies and find the dual
                             !! Egrid different from zero, with psp_modifier=0,  activates the fit of psigrid
                             !! energy by energy and the calculation of paw stuff
  print *, "copy psigrid_pseudo "
  psigrid_pseudo=psigrid
  !!$  print *, "  chiamo   abs_generator_modified con psipsigrid_pseudo, Ngrid_biggerbox  ", Ngrid_biggerbox
  !!$ psipsigrid_pseudo(1,1)=Ngrid_biggerbox
  call abs_generator_modified(iproc,atoms%nzatom(ity), atoms%nelpsp(ity),atoms%psppar(:,:,ity),&
       atoms%npspcode(ity),ng-1 ,noccmax , lmax , expo,psi,aeval, occup , psp_modifier , &
       Nsol, abs_final_L , Ngrid, Ngrid_box,Egrid_pseudo,  rgrid , psigrid_pseudo ,&
       Npaw, PAWpatch,  psipsigrid_pseudo)
  
  if(dump_functions==1)  then
     if(iproc.eq.0) then
        do iocc=1,2
           print *, " pseudo  Egau, pseudo  Egrid ", aeval(iocc, abs_final_L+1  ), Egrid_pseudo(iocc)
        enddo
     endif
  endif
  if(iproc.eq.0)   print * , "  routine GetExcitedOrbitalAsG  , calculate pseudo OK  " 
  if(iproc.eq.0 .and. dump_functions.eq.1) then
     if(psp_modifier.eq.0) then
        call dump_gauwf_on_radgrid("pseudo_wf_radgrid", ng-1,noccmax,lmax,expo,psi)
     else
        call dump_gauwf_on_radgrid("real_wf_radgrid",   ng-1,noccmax,lmax,expo,psi)
     endif
  endif
  if(iproc.eq.0 .and. dump_functions.eq.1) then
     open(unit=22,file='numerov_ptildes.dat')
     do igrid=1, Ngrid
        write(22,'(200(f20.10,1x))') rgrid(igrid), (psigrid_pseudo(igrid,j ), j=1,Npaw)   
     enddo
     close(unit=22)
     open(unit=22,file='numerov_psitildes.dat')
     do igrid=1, Ngrid
        write(22,'(200(f20.10,1x))') rgrid(igrid), (psipsigrid_pseudo(igrid,j ), j=1,Npaw)   
     enddo
     close(unit=22)
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
  

  if(.true. .or. dump_functions==1)  then
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
  endif
  

  !!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !!  In the following   AE are calculated on a bigger grid to 
  !! avoid compression of the  missing bound states
  !!
  if(iproc.eq.0) print * , " routine GetExcitedOrbitalAsG  , calculate  atom with  real-pot, abs_final_L " ,abs_final_L
  psp_modifier=1
  call abs_generator_modified(iproc,atoms%nzatom(ity), atoms%nelpsp(ity),atoms%psppar(:,:,ity),&
       atoms%npspcode(ity),Ngrid ,noccmax , lmax , expo,psireal,aeval, occup , psp_modifier ,  &
       Nsol, abs_final_L , Ngrid, Ngrid_biggerbox,Egrid,  rgrid , psigrid_bigger ,Npaw,  PAWpatch )

  if(iproc.eq.0 .and. dump_functions.eq.1) then
     open(unit=22,file='numerov_bigger.dat')
     do igrid=1, Ngrid
        write(22,'(200(f20.10,1x))') rgrid(igrid), (psigrid_bigger(igrid,j ), j=1,Nsol)   
     enddo
     close(unit=22)
  endif



  if(iproc.eq.0) print *, " routine GetExcitedOrbitalAsG  , PROJECT 1s*r on pseudos "
  print *, "   shape(coeffs_in_paw) "  ,  shape(coeffs_in_paw)
  call find_pfproj_4tail( Nsol,Npaw,Ngrid,  Ngrid_box,Ngrid_biggerbox, rgrid, psi1s, psigrid, real_start, &
       psigrid_pseudo, psipsigrid_pseudo,  &
       psigrid_bigger,dump_functions, coeffs_in_paw)


  !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!!$  call find_pfproj( Nsol,Nsol-5,Ngrid,  Ngrid_biggerbox, rgrid, psi1s, psigrid, real_start, &
!!$       psipsigrid_pseudo, dump_functions, coeffs_in_paw)

  write(filename_in_coeffs,'(A,A)') trim(filename), "_patches"
  open(unit=22,file=filename_in_coeffs, status='unknown')
     print *, "now writing the coefficients of the initial wave function in psitilde basis  "
     write(22,'(A)') " now writing the coefficients of the initial wave function in psitilde basis "
 
     write(22,'(I4)') Npaw
     do isol=1, Npaw
        write(22,'(E20.14)') coeffs_in_paw(isol)   
     end do
  close(unit=22)

  print * , "  OK pfproj " 

  !! if(iproc.eq.0 .and. dump_functions.eq.1) then 
  if(iproc.eq.0 ) then 
     open(unit=22,file='projres.dat')
     do igrid=1, Ngrid_biggerbox
        !! la funzione proiettata e in ultima colonna
        write(22,'(200(f20.10,1x))') rgrid(igrid),  psi1s(igrid), psigrid(igrid,1), psigrid(igrid,2)
     enddo

!!$     sum=0.0
!!$     do igrid=1, Ngrid-1
!!$        sum=sum+(rgrid(igrid+1)-rgrid(igrid)) *psigrid(igrid,2) **2
!!$     end do
!!$     print *, "  norma 2 del proiettato ", sum

     close(unit=22)
  endif
  


!!$  expo_fine(1)=             expo(1)/15.0
!!$  expo_fine(ng_fine) =       cradius/2
!!$  
!!$  do ig=1, ng_fine
!!$     expo_fine(ig) = exp( ( log(expo_fine(1))*(ng_fine-ig) +log(expo_fine(ng_fine))*(ig-1) )/(ng_fine-1))
!!$  enddo
!!$  
!!$  call find_Scoeffs_grid(ng_fine-1,  expo_fine, Ngrid, rgrid, psigrid(:,2)  , gcoeffs , abs_final_L  )
!!$
!!$  !! if(iproc.eq.0 .and. dump_functions.eq.1)  call  dump_1gauwf_on_radgrid("proje_gau_proje_pseudo.dat",&
!!$  if(iproc.eq.0 )  call  dump_1gauwf_on_radgrid("proje_gau_proje_pseudo.dat",&
!!$       ng_fine-1 , expo_fine,gcoeffs   ,   abs_final_L +1  )
!!$  
!!$ 
!!$  
!!$  ! ----------------- Gabsorber --------------------------------------------------------
!!$  Gabsorber%nat = 1
!!$  allocate(Gabsorber%rxyz(3,Gabsorber%nat+ndebug),stat=i_stat)
!!$  call memocc(i_stat,Gabsorber%rxyz ,'Gabsorber%rxyz',subname)
!!$
!!$  allocate(Gabsorber%nshell(Gabsorber%nat+ndebug),stat=i_stat)
!!$  call memocc(i_stat,Gabsorber%nshell,'Gabsorber%nshell',subname)
!!$
!!$  Gabsorber%rxyz(:,1) = rxyz(:,in_iat_absorber )
!!$
!!$  Gabsorber%nshell(1)=1
!!$  Gabsorber%nshltot  =1
!!$
!!$  allocate(Gabsorber%ndoc(1+ndebug),stat=i_stat)
!!$  call memocc(i_stat,   Gabsorber%ndoc,    'Gabsorber%ndoc',subname)
!!$  allocate(Gabsorber%nam (1+ndebug),stat=i_stat)
!!$  call memocc(i_stat,Gabsorber%nam,'Gabsorber%nam',subname)
!!$
!!$  Gabsorber%nexpo=0
!!$  Gabsorber%ncoeff=0
!!$
!!$  Gabsorber%ndoc (1)  =  ng_fine
!!$  Gabsorber%nam  (1)   =  abs_final_l+1
!!$  Gabsorber%nexpo         =  Gabsorber%nexpo+ ng_fine 
!!$  Gabsorber%ncoeff        =  Gabsorber%ncoeff+2* abs_final_l + 1
!!$
!!$  allocate(Gabsorber%psiat(Gabsorber%nexpo+ndebug),stat=i_stat)
!!$  call memocc(i_stat,Gabsorber%psiat , 'Gabsorber%psiat',subname)
!!$
!!$  allocate(Gabsorber%xp(Gabsorber%nexpo+ndebug),stat=i_stat)
!!$  call memocc(i_stat,Gabsorber%xp    , 'Gabsorber%xp',subname)
!!$
!!$  iexpo=0
!!$  do ig=1,Gabsorber%ndoc(1)
!!$     iexpo=iexpo+1
!!$     Gabsorber%psiat(iexpo)=gcoeffs (ig-1)
!!$     Gabsorber%xp(iexpo)=expo_fine(ig)
!!$  end do
!!$
!!$  print *,"writing gaussian form of the projected initial wave on file " , filename
!!$  !! psigrid is used just to fill the argument box.
!!$  !! The coefficient are alway reinitialised according to polarisation
!!$  call write_gaussian_information( 0 ,1 ,ha%orbs,Gabsorber,psigrid  ,filename)


  if(in_NPaw.gt.0 ) then
  
     open(unit=22,file=filename_in_coeffs, access='append')
     
     print *, "psi1s : writing the projected initial wf "
       write(22,'(A,A)') "psi1s : writing the projected initial wf on file  " , trim(filename_in_coeffs)
       write(22,'(I7)')  Ngrid_biggerbox
       do igrid=1, Ngrid_biggerbox 
          !! write(22,'(E20.14,1x,E20.14)') rgrid(igrid), psi1s(igrid)
          write(22,'(E20.14,1x,E20.14)') rgrid(igrid),  psigrid(igrid,2)
       enddo

     print *, "now writing the dual functions and PAWpatch for ", Npaw," duals "
       write(22,'(A)') "now writing the dual functions and PAWpatch  "
       write(22,'(I4,1x,I6,1x,I4)') Npaw, Ngrid_box, abs_final_L
       do igrid=1, Ngrid_box
          write(22,'(E20.14)') rgrid(igrid)
       enddo
     
     write(22,'(A)') " "
     do isol=1, Npaw
        do igrid=1, Ngrid_box
           write(22,'(E20.14)') psigrid_pseudo(igrid, isol)
        enddo
        write(22,'(A)') " "
     enddo
     
     do isol=1, Npaw
        do jsol=1, Npaw
           write(22,'(E20.14)') PawPatch(jsol,isol)
        enddo
     enddo
     
     close(unit=22)

  endif
  





  write(filename_in_coeffs,'(A,A)') trim(filename), "_initial"
  open(unit=22,file=filename_in_coeffs, status='unknown')
     print *, "now writing the coefficients of the initial wave function in psitilde basis  "
     write(22,'(A)') " now writing the coefficients of the initial wave function in psitilde basis "
 
     write(22,'(I4)') 1
     do isol=1, 1
        write(22,'(E20.14)') 1.0_8   
     end do
  close(unit=22)
  


  if(in_NPaw.gt.0 ) then
  
     open(unit=22,file=filename_in_coeffs, access='append')
 
     print *, "now writing the dual functions and PAWpatch for ", 1," duals "
     write(22,'(A)') "now writing the dual functions and PAWpatch  "
     write(22,'(I4,1x,I6,1x,I4)') 1, Ngrid_biggerbox, -abs_final_L-1
     do igrid=1, Ngrid_biggerbox
        write(22,'(E20.14)') rgrid(igrid)
     enddo
     
     write(22,'(A)') " "
     do isol=1, 1
        do igrid=1, Ngrid_biggerbox
           write(22,'(E20.14)') psigrid(igrid,2)
        enddo
        write(22,'(A)') " "
     enddo
     
     do isol=1, 1
        do jsol=1, 1
           write(22,'(E20.14)') 0.0_8
        enddo
     enddo
     
     close(unit=22)

  endif
  





















  ! print *,'expo',shape(expo),ng_fine,expo(:)

  ! -------------------------------------------------------------------------
!!$
!!$


  i_all=-product(shape(PAWpatch))*kind(PAWpatch)
  deallocate(PAWpatch,stat=i_stat)
  call memocc(i_stat,i_all,'PAWpatch',subname)



  i_all=-product(shape(psigrid_pseudo))*kind(psigrid_pseudo)
  deallocate(psigrid_pseudo,stat=i_stat)
  call memocc(i_stat,i_all,'psigrid_pseudo',subname)
  


  i_all=-product(shape(psipsigrid_pseudo))*kind(psipsigrid_pseudo)
  deallocate(psipsigrid_pseudo,stat=i_stat)
  call memocc(i_stat,i_all,'psipsigrid_pseudo',subname)
  

  i_all=-product(shape(psigrid_bigger))*kind(psigrid_bigger)
  deallocate(psigrid_bigger,stat=i_stat)
  call memocc(i_stat,i_all,'psigrid_bigger',subname)
  


  i_all=-product(shape(psigrid))*kind(psigrid)
  deallocate(psigrid,stat=i_stat)
  call memocc(i_stat,i_all,'psigrid',subname)
  
  i_all=-product(shape(rgrid))*kind(rgrid)
  deallocate(rgrid,stat=i_stat)
  call memocc(i_stat,i_all,'rgrid',subname)




  i_all=-product(shape(Egrid))*kind(Egrid)
  deallocate(Egrid,stat=i_stat)
  call memocc(i_stat,i_all,'Egrid',subname)


 
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

  i_all=-product(shape(psireal))*kind(psireal)
  deallocate(psireal,stat=i_stat)
  call memocc(i_stat,i_all,'psireal',subname) 
  
  i_all=-product(shape(expo))*kind(expo)
  deallocate(expo,stat=i_stat)
  call memocc(i_stat,i_all,'expo',subname)
    
  i_all=-product(shape(expo_fine))*kind(expo_fine)
  deallocate(expo_fine,stat=i_stat)
  call memocc(i_stat,i_all,'expo_fine',subname)
  
  i_all=-product(shape(Egrid_pseudo))*kind(Egrid_pseudo)
  deallocate(Egrid_pseudo,stat=i_stat)
  call memocc(i_stat,i_all,'Egrid_pseudo',subname)

  i_all=-product(shape(psi1s))*kind(psi1s)
  deallocate(psi1s,stat=i_stat)
  call memocc(i_stat,i_all,'psi1s',subname)
  
  if(dump_functions==1)  then
     !! stop
  end if

  return
END SUBROUTINE GetExcitedOrbitalAsG


function GetBottom( atoms, iproc)
  
  use module_base
  use module_types
  use module_interfaces,except_this_one => GetBottom

  implicit none
  !Arguments
  type(atoms_data), intent(in) :: atoms
  real(gp) :: GetBottom
  integer :: iproc
  !Local variables
  character(len=*), parameter :: subname='GetBottom'
  integer :: abs_final_L
  integer :: ity, ng , noccmax, lmax, i_all
  real(gp) , pointer :: expo(:), psi(:,:,:), aeval(:,:), occup(:,:)
  integer :: psp_modifier

  integer, parameter :: Norder=4, dump_functions=0
  real(gp) :: cradius

  integer ::  Nsol , Ngrid, igrid

  real(gp), pointer :: Egrid(:) ,  rgrid(:) , psigrid (:,:)
  integer :: i_stat, Npaw
  real(gp) :: rzero
  real(gp)  dumPAWpatch(1,1)  !! not used here. Just to fill an argument

  ! if (in_iat_absorber.ne.0) then

  ng  = 21
  noccmax = 5 
  lmax=3

  Nsol=2
  Ngrid=3000

  cradius=5.0 !!!!!!!! ATTENZIONE

  allocate(expo(ng +ndebug  ), stat=i_stat)
  call memocc(i_stat,expo,'expo',subname)

  allocate(psi ( 0:ng-1  ,noccmax,lmax+1+ndebug ), stat=i_stat)
  call memocc(i_stat,psi,'psi',subname)

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
  

  rzero = 1.0D-2/Ngrid * cradius 
  do igrid=1, Ngrid
     rgrid(igrid) = rzero*  exp( igrid*   1.0_gp/Ngrid * log( cradius/rzero ))
  enddo
  
  abs_final_L = 0

  GetBottom=1.0D4

  psp_modifier=0;

  do ity=1, atoms%ntypes
     Egrid=0.0_gp
     !! to fit these energies and find the dual
     !! Egrid different from zero activates the fit of psigrid
     !! energy by energy and the calculation of paw stuff
     Npaw=1
     call abs_generator_modified(iproc,atoms%nzatom(ity), atoms%nelpsp(ity),atoms%psppar(:,:,ity),&
          atoms%npspcode(ity),ng-1 ,noccmax , lmax , expo,psi,aeval, occup , psp_modifier , &
          Nsol, abs_final_L , Ngrid,Ngrid, Egrid,  rgrid , psigrid, Npaw, dumPAWpatch )
     ! if(aeval(1,1)<GetBottom) GetBottom=aeval(1,1)
     !! print *,  atoms%atomnames(ity)
     !! print *, aeval
     if( minval(aeval) <GetBottom) GetBottom=minval(aeval)
  enddo

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

  i_all=-product(shape(psi))*kind(psi)
  deallocate(psi,stat=i_stat)
  call memocc(i_stat,i_all,'psi',subname) 

  i_all=-product(shape(expo))*kind(expo)
  deallocate(expo,stat=i_stat)
  call memocc(i_stat,i_all,'expo',subname)

  i_all=-product(shape(rgrid))*kind(rgrid)
  deallocate(rgrid,stat=i_stat)
  call memocc(i_stat,i_all,'rgrid',subname)

end function GetBottom


subroutine zero4b2B(n,x)
  implicit none
  !Arguments
  integer, intent(in) :: n
  real(kind=8), intent(out) :: x(n)
  !Local variables
  integer :: i
  do i=1,n
     x(i)=0.d0
  end do
END SUBROUTINE zero4b2B


!> backward wavelet transform
!! nd: length of data set
!! nt length of data in data set to be transformed
!! m filter length (m has to be even!)
!! x input data, y output data
subroutine back_trans_14_4b2B(nd,nt,x,y)
  implicit none
  !Arguments
  integer, intent(in) :: nd,nt
  real(kind=8), intent(in) :: x(0:nd-1)
  real(kind=8), intent(out) :: y(0:nd-1)
  !Local variables
  integer :: i,j,ind

  include 'lazy_16.inc'
  
  do i=0,nt/2-1
     y(2*i+0)=0.d0
     y(2*i+1)=0.d0
     
     do j=-m/2,m/2-1
        
        ! periodically wrap index if necessary
        ind=i-j
        loop99: do
           if (ind.lt.0) then 
              ind=ind+nt/2
              cycle loop99
           end if
           if (ind.ge.nt/2) then 
              ind=ind-nt/2
              cycle loop99
           end if
           exit loop99
        end do loop99

        y(2*i+0)=y(2*i+0) + ch(2*j-0)*x(ind)+cg(2*j-0)*x(ind+nt/2)
        y(2*i+1)=y(2*i+1) + ch(2*j+1)*x(ind)+cg(2*j+1)*x(ind+nt/2)
     end do
  end do

END SUBROUTINE back_trans_14_4b2B


subroutine scaling_function4b2B(itype,nd,nrange,a,x)
  use module_base
  implicit none
  !Arguments
  !Type of interpolating functions
  integer, intent(in) :: itype
  !Number of points: must be 2**nex
  integer, intent(in) :: nd
  integer, intent(out) :: nrange
  real(kind=8), dimension(0:nd), intent(out) :: a,x
  !Local variables
  character(len=*), parameter :: subname='scaling_function4b2B'
  real(kind=8), dimension(:), allocatable :: y
  integer :: i,nt,ni,i_all,i_stat  

  !Only itype=8,14,16,20,24,30,40,50,60,100
  select case(itype)
  case(8,14,16,20,24,30,40,50,60,100)
     !O.K.
  case default
     print *,"Only interpolating functions 8, 14, 16, 20, 24, 30, 40, 50, 60, 100"
     stop
  end select
!!$  write(unit=*,fmt="(1x,a,i0,a)") &
!!$       "Use interpolating scaling functions of ",itype," order"

  !Give the range of the scaling function
  !from -itype to itype
  ni=2*itype
  nrange = ni
  allocate(y(0:nd+ndebug),stat=i_stat)
  call memocc(i_stat,y,'y',subname)
  
  ! plot scaling function
  call zero4b2B(nd+1,x)
  call zero4b2B(nd+1,y)
  nt=ni
  x(nt/2)=1.d0
  loop1: do
     nt=2*nt
     ! write(6,*) 'nd,nt',nd,nt
     select case(itype)
     case(8)
        stop
     case(14)
        stop
     case(16)
        call back_trans_14_4b2B(nd,nt,x,y)
     case(20)
        stop
     case(24)
        stop
     case(30)
        stop
     case(40)
        stop
     case(50)
        stop
     case(60)
        stop
     case(100)
        stop
     end select

     do i=0,nt-1
        x(i)=y(i)
     end do
     if (nt.eq.nd) then
        exit loop1
     end if
  end do loop1

  !open (unit=1,file='scfunction',status='unknown')
  do i=0,nd
     a(i) = real(i*ni,kind=8)/real(nd,kind=8)-(.5d0*real(ni,kind=8)-1.d0)
     !write(1,*) a(i),x(i)
  end do
  !close(1)

  i_all=-product(shape(y))*kind(y)
  deallocate(y,stat=i_stat)
  call memocc(i_stat,i_all,'y',subname)
END SUBROUTINE scaling_function4b2B


subroutine read_potfile4b2B(filename,n1i,n2i,n3i, rho, alat1, alat2, alat3)
  use module_base
  implicit none
  character(len=*), intent(in) :: filename
  integer, intent(out) :: n1i,n2i,n3i
  real(gp) alat1, alat2, alat3, dum, dum1
  ! real(dp), dimension(n1i*n2i*n3d), intent(out) :: rho
  real(gp), pointer :: rho(:)
  !local variables
  integer :: nl1,nl2,nl3,i_stat,i1,i2,i3,ind
  real(gp) :: value
  character(len=*), parameter :: subname='read_potfile4b2B'

  open(unit=22,file=filename,status='unknown')
  read(22,*)!'normalised density'
  read(22,*) n1i,n2i,n3i
  read(22,*) alat1,dum ,alat2
  read(22,*)  dum, dum1, alat3
  read(22,*)!xyz   periodic' !not true in general but needed in the case

  !conditions for periodicity in the three directions
  !value of the buffer in the x and z direction
  nl1=1
  nl3=1
  nl2=1

  print *, " allocation for rho for  n1i,n2i,n3i ",  n1i,n2i,n3i

  allocate( rho( n1i*n2i*n3i+ndebug) , stat=i_stat )
  call memocc(i_stat,rho,'rho',subname)

  print *, " going to read all pot points " 
  do i3=0,n3i-1
     do i2=0,n2i-1
        do i1=0,n1i-1
           ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
           read(22,*)value
           rho(ind)=value
        end do
     end do
  end do
  print *, " closing file  " 
  close(22)
  
END SUBROUTINE read_potfile4b2B



subroutine  atom_numeric(rprb,lmax,noccmax,occup,&
                 zion,&
                 aeval,psi,&
                 Nsol, Labs, Ngrid,Ngrid_box, Egrid,  rgrid , psigrid )
  use module_base, only: gp
  use esatto
  implicit none
  
  
  real(gp) :: psi(Ngrid ,noccmax,lmax+1),aeval(noccmax,lmax+1),&
       occup(noccmax,lmax+1),&
       psigrid(Ngrid, Nsol), &
       rhogrid(Ngrid), potgrid(Ngrid), rhogrid_old(Ngrid), &
       vxcgrid(Ngrid), &
       Egrid(nsol)
  
  real(gp) :: rgrid(Ngrid), rprb, zion
  real(gp), target :: dumgrid1(Ngrid),dumgrid2(Ngrid)

  integer :: igrid, l, iocc, i, it, isol, Labs, Ngrid, lmax, Nsol, noccmax, Ngrid_box
  real(gp) :: evsumold, evsum, rmix, tt, dum, r
  real(gp) :: emuxc, totocc


  print *, "INSIDE  ATOM_NUMERIC , RPRB ", rprb, " Nsol ", Nsol


  do igrid=1,Ngrid
     rhogrid(igrid)=0.0_gp
  enddo

  evsum=1.d30

  big_loop: do it=1,200
     evsumold=evsum
     evsum=0._gp
     
! coefficients of charge density
     rhogrid_old(:)=rhogrid(:)
     rhogrid(:)=0.0_gp
     totocc=0.0_gp
     do l=0,lmax
        do iocc=1,noccmax
           if (occup(iocc,l+1).gt.0._gp) then
              totocc=totocc+occup(iocc,l+1)
              do i=1,Ngrid
                    rhogrid(i)=rhogrid(i) + &
                         psi(i,iocc,l+1)*psi(i,iocc,l+1)*occup(iocc,l+1)
              end do
           end if
        end do
     end do



     rmix=.5_gp
     if (it.eq.1) rmix=1._gp
     do i=1,Ngrid
        tt=rmix*rhogrid(i) + (1._gp-rmix)*rhogrid_old(i)
        rhogrid(i)=tt
     end do
     rhogrid_old(:)=rhogrid(:)



     do igrid=1,Ngrid
        r=rgrid(igrid)
        dum = rhogrid(igrid)/r/r *0.07957747154594768_gp
        vxcgrid(igrid)=emuxc(dum) 
     enddo
     ! print *, rgrid(:10), potgrid(:10)

     do igrid=1, Ngrid
        r=rgrid(igrid)
        potgrid(igrid) =0.5_gp*r*r  /    rprb**4 
        potgrid(igrid) = potgrid(igrid) - zion/r
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
        do igrid=1, ngrid
           r=rgrid(igrid)
           dumgrid2(igrid)=potgrid(igrid)+ 0.5_gp*l*(l+1.0_gp)/r/r + vxcgrid(igrid)
        enddo
        
        do iocc=1,noccmax
           if (occup(iocc,l+1).gt.0._gp) then
              
              
              dumgrid1(:)=0.0_gp
              call schro(aeval(iocc,l+1),rgrid(1),dumgrid2(1),dumgrid1(1), psi(1,iocc,l+1) ,ngrid, iocc+l,l*1.0_gp, zion)
              evsum=evsum+aeval(iocc,l+1)

             
           end if
        end do
     end do loop_l
        
     tt=abs(evsum-evsumold)

     if (tt.lt.1.e-12_gp) then
           exit big_loop
     end if
  end do big_loop
  ! End of the big loop


  do igrid=1, ngrid
     r=rgrid(igrid)
     potgrid(igrid)=potgrid(igrid)+ 0.5_gp*labs*(labs+1.0_gp)/r/r + vxcgrid(igrid)
  enddo
  
  dumgrid1(:)=0.0_gp
  do isol=1,nsol
     psigrid(:,isol)=0.0_gp
     call schro(Egrid(isol),rgrid(1),potgrid(1),dumgrid1(1),psigrid(1,isol),ngrid_box,isol+labs,labs*1.0_gp, zion)

  enddo

END SUBROUTINE atom_numeric





subroutine  atom_numeric_eqdiff(rprb,lmax,noccmax,occup,&
                 zion,&
                 aeval,psi,&
                 Nsol, Labs, Ngrid,Egrid,  rgrid , psigrid )
  use module_base
  use esatto
  implicit none
  
  
  real(gp) :: psi(Ngrid ,noccmax,lmax+1),aeval(noccmax,lmax+1),&
       occup(noccmax,lmax+1),&
       psigrid(Ngrid, Nsol), &
       rhogrid(Ngrid), potgrid(Ngrid), rhogrid_old(Ngrid), &
       vxcgrid(Ngrid), &
       Egrid(nsol)
  


  real(gp) ::  zion, phase, fattore
  real(gp) rprb
  integer :: igrid, l, iocc, i, it, isol, Labs, Ngrid, lmax, Nsol, noccmax, idiff, Nrdiff, nls_a, i_stat
  real(gp) :: evsumold, evsum, rmix, tt, dum, r
  real(gp) :: emuxc, totocc
  real(gp) Ediff, Rdiff, rpot_a,spot_a,hpot_a, Rinf_a

  real(gp) :: rgrid(Ngrid)
  real(gp), target :: dumgrid1(Ngrid),dumgrid2(Ngrid) ,dumgrid3(Ngrid)
  real(gp) , pointer :: y_r(:) ,  d_r(:)
  real(gp) :: ref
  integer :: nsteps_fine, nstesp_coarse, lpot_a, i_all
  character(len=*), parameter :: subname='atom_numeric_eqdiff'



  do igrid=1,Ngrid
     rhogrid(igrid)=0.0_gp
  enddo

  evsum=1.d30

  big_loop: do it=1,200
     evsumold=evsum
     evsum=0._gp
     
! coefficients of charge density
     rhogrid_old(:)=rhogrid(:)
     rhogrid(:)=0.0_gp
     totocc=0.0_gp
     do l=0,lmax
        do iocc=1,noccmax
           if (occup(iocc,l+1).gt.0._gp) then
              totocc=totocc+occup(iocc,l+1)
              do i=1,Ngrid
                    rhogrid(i)=rhogrid(i) + &
                         psi(i,iocc,l+1)*psi(i,iocc,l+1)*occup(iocc,l+1)
              end do
           end if
        end do
     end do



     rmix=.5_gp
     if (it.eq.1) rmix=1._gp
     do i=1,Ngrid
        tt=rmix*rhogrid(i) + (1._gp-rmix)*rhogrid_old(i)
        rhogrid(i)=tt
     end do
     rhogrid_old(:)=rhogrid(:)



     do igrid=1,Ngrid
        r=rgrid(igrid)
        dum = rhogrid(igrid)/r/r *0.07957747154594768_gp
        vxcgrid(igrid)=emuxc(dum) 
     enddo
     ! print *, rgrid(:10), potgrid(:10)

     do igrid=1, Ngrid
        r=rgrid(igrid)
        potgrid(igrid) =0.5_gp*r*r  /    rprb**4 
        potgrid(igrid) = potgrid(igrid) - zion/r
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
        do igrid=1, ngrid
           r=rgrid(igrid)
           dumgrid2(igrid)=potgrid(igrid)+ 0.5_gp*l*(l+1.0_gp)/r/r + vxcgrid(igrid)
        enddo
        
        do iocc=1,noccmax
           if (occup(iocc,l+1).gt.0._gp) then
              
              
              dumgrid1(:)=0.0_gp
              call schro(aeval(iocc,l+1),rgrid(1),dumgrid2(1),dumgrid1(1), psi(1,iocc,l+1) ,ngrid, iocc+l,l*1.0_gp, zion)
              evsum=evsum+aeval(iocc,l+1)

             
           end if
        end do
     end do loop_l

     print *," evsum " ,  evsum
     tt=abs(evsum-evsumold)

     if (tt.lt.1.e-12_gp) then
           exit big_loop
     end if
  end do big_loop

  ! End of the big loop
  do igrid=1, Ngrid
     r=rgrid(igrid)
     potgrid(igrid)=potgrid(igrid) - 0.5_gp*r*r  /    rprb**4 
  enddo

  open(unit=22,file='pot.dat')
  do igrid=1, Ngrid
     write(22,'(200(E20.10,1x))') rgrid(igrid),potgrid(igrid)+vxcgrid(igrid), potgrid(igrid)
  enddo
  close(unit=22)

  nls_a=16

  allocate(y_r(0:nls_a-1+ndebug), stat=i_stat)
  call memocc(i_stat,y_r,'y_r',subname)
  allocate(d_r(0:nls_a-1+ndebug), stat=i_stat)
  call memocc(i_stat,d_r,'d_r',subname)
  

  
  
     
  print *, " trying  eq diff "
  do idiff=1,600
     Ediff=idiff*0.04_gp/3.0 +0.0_gp
     do igrid=Ngrid,1,-1
        dumgrid2(igrid)=potgrid(igrid)+vxcgrid(igrid)
        if(rgrid(igrid)>3.5_gp) then
           Nrdiff = igrid;
           Rdiff = rgrid(igrid)
        endif
     enddo
     do l=0, nls_a-1
        if(l>0) then
           do igrid=1,Ngrid
              dumgrid2 (igrid)=dumgrid2 (igrid)+ 0.5_gp*(2.0_gp*l)/rgrid(igrid)/rgrid(igrid)
           enddo
        endif
        dumgrid1=0.0_wp
        d_r(l)=phase( Ediff, Nrdiff, rgrid,dumgrid2  , dumgrid1 , psigrid(1,2+l) , l*1.0_gp ,0, 1)
        y_r(l)= psigrid(Nrdiff,2+l)
        
     enddo
     
     dumgrid1 = 0.0_gp
     dumgrid2 = potgrid + vxcgrid
     call schro(Egrid(1) , rgrid , dumgrid2  , dumgrid1, psigrid(:,1) , ngrid , 1+0 , 0*1.0_gp ,  zion)
     
        
     open(unit=22,file='bound.dat')
     do igrid=1, Ngrid
        write(22,'(200(f20.10,1x))') rgrid(igrid),  psigrid(igrid,1), psigrid(igrid,2+Labs)
     enddo
     close(unit=22)
     

     dumgrid1(:) =  psigrid(:,2+Labs)* psigrid( :, 1)
     do l=1,Labs
        dumgrid1(:) = dumgrid1(:) *rgrid(:)
     end do
     call integrate(dumgrid1(1),dumgrid3(1),rgrid(1) ,Nrdiff)
     fattore = dumgrid3(Nrdiff)

     dumgrid1(:)=0
     if(.true.) then
        lpot_a=1
        lpot_a=1
        rpot_a = 7.5_gp
        spot_a = 0.8_gp
        hpot_a = 3.0_gp
        Rinf_a=100.0_gp
        nstesp_coarse=1000
        nsteps_fine  = 40
        ref= esatto_CalcolaRiflettivita( ngrid, rgrid, dumgrid2, nls_a, lpot_a, rpot_a,spot_a,hpot_a,y_r,d_r,&
             Rdiff,    Rinf_a ,nstesp_coarse ,nsteps_fine, Ediff , Labs)
        
        print *," ",  Ediff, ref*fattore*fattore, fattore
     endif
  enddo
  i_all=-product(shape(y_r))*kind(y_r)
  deallocate(y_r, stat=i_stat)
  call memocc(i_stat,i_all,'y_r',subname)
  i_all=-product(shape(d_r))*kind(d_r)
  deallocate(d_r, stat=i_stat)
  call memocc(i_stat,i_all,'d_r',subname)
  stop
END SUBROUTINE atom_numeric_eqdiff
