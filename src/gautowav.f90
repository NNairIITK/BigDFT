subroutine gautowav(iproc,nproc,nat,ntypes,norb,norbp,n1,n2,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,&
     nvctr_c,nvctr_f,nseg_c,nseg_f,keyg,keyv,iatype,occup,rxyz,hgrid,psi,eks)

  implicit none
  integer, intent(in) :: norb,norbp,iproc,nproc,nat,ntypes
  integer, intent(in) :: nvctr_c,nvctr_f,n1,n2,n3,nseg_c,nseg_f
  integer, intent(in) :: nfl1,nfu1,nfl2,nfu2,nfl3,nfu3
  integer, dimension(nseg_c+nseg_f), intent(in) :: keyv
  integer, dimension(2,nseg_c+nseg_f), intent(in) :: keyg
  integer, dimension(nat), intent(in) :: iatype
  real(kind=8), intent(in) :: hgrid
  real(kind=8), intent(in) :: rxyz(3,nat)
  real(kind=8), dimension(norb), intent(in) :: occup
  real(kind=8), intent(out) :: eks
  real(kind=8), dimension(nvctr_c+7*nvctr_f,norbp), intent(out) :: psi
  !local variables
  logical :: myorbital
  character(len=6) :: string,symbol
  character(len=100) :: line
  integer, parameter :: nterm_max=3
  integer :: ngx,nbx,npgf,nst,nend,ng,lshell,num,mmax,myshift,icbas,isbas,nbas,nco,i,ipar,ipg,jat
  integer :: iorb,jorb,iat,ityp,l,m,nterm,i_all,i_stat,ibas,ig,iset,jbas,iterm,ishell,lmax
  real(kind=8) :: rx,ry,rz,anorm,coeff,const0,const1,exponent,coefficient,scpr,ek,tt
  integer, dimension(:), allocatable :: lx,ly,lz,nshell,iorbtmp
  integer, dimension(:,:), allocatable :: nam,ndoc
  real(kind=8), dimension(:), allocatable :: fac_arr,psiatn,xp,tpsi,ctmp
  real(kind=8), dimension(:,:,:), allocatable :: contcoeff,expo
  real(kind=8), dimension(:,:,:,:), allocatable :: cimu

  !parse the output of CP2K to read the basis set information

  ngx=0
  nbx=0
  lmax=0

  allocate(nshell(ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(nshell))*kind(nshell),'nshell','gautowav')

  open(unit=35,file='gaubasis.dat',action='read')

  !read the lines for analyzing the maximum number of primitive gaussian functions
  !and als for the maximum angular momentum
  ityp=0
  ngx=0
  nbx=0
  lmax=0
  ipg=0
  for_ngx: do
     if (ityp > ntypes) exit for_ngx
     read(35,'(a100)')line
     !analyzing the different possibilities
     read(line,*,iostat=i_stat)tt,string,symbol
     if (i_stat == 0 .and. string=='Atomic' .and. symbol=='kind:') then
        ityp=ityp+1
        if (ityp > 1) then
           nshell(ityp-1)=ishell
           nbx=max(nbx,ishell)
        end if
        ishell=0
        cycle for_ngx
     end if
     read(line,*,iostat=i_stat)iset,num,num,num,exponent,coefficient
     if (i_stat==0) then
        !print *,num,exponent,coefficient
        ishell=ishell+1
        lmax=max(lmax,num)
        ngx=max(ngx,ipg)
        !print *,ishell,ipg,lmax
        ipg=1
        cycle for_ngx
     end if
     read(line,*,iostat=i_stat)exponent,coefficient
     if (i_stat==0) then
        ipg=ipg+1
        cycle for_ngx
     end if
  end do for_ngx

  !now store the values
  rewind(35)

  !here allocate arrays
  allocate(nam(nbx,ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(nam))*kind(nam),'nam','gautowav')
  allocate(ndoc(nbx,ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(ndoc))*kind(ndoc),'ndoc','gautowav')
  allocate(contcoeff(ngx,nbx,ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(contcoeff))*kind(contcoeff),'contcoeff','gautowav')
  allocate(expo(ngx,nbx,ntypes),stat=i_stat)
  call memocc(i_stat,product(shape(expo))*kind(expo),'expo','gautowav')
  
  ityp=0
  ipg=0
  store_basis: do
     if (ityp > ntypes) exit store_basis
     read(35,'(a100)')line
     !analyzing the different possibilities
     read(line,*,iostat=i_stat)tt,string,symbol
     if (i_stat == 0 .and. string=='Atomic' .and. symbol=='kind:') then
        ityp=ityp+1
        if (ityp > 1) then 
           ndoc(ishell,ityp-1)=ipg
        end if
        ishell=0
        cycle store_basis
     end if
     read(line,*,iostat=i_stat)iset,num,num,num,exponent,coefficient
     if (i_stat==0) then
        !print *,num,exponent,coefficient
        ishell=ishell+1
        nam(ishell,ityp)=num
        lmax=max(lmax,num)
        if (ishell > 1) ndoc(ishell-1,ityp)=ipg
        expo(1,ishell,ityp)=exponent
        contcoeff(1,ishell,ityp)=coefficient
        ipg=1
        cycle store_basis
     end if
     read(line,*,iostat=i_stat)exponent,coefficient
     if (i_stat==0) then
        ipg=ipg+1
        expo(ipg,ishell,ityp)=exponent
        contcoeff(ipg,ishell,ityp)=coefficient
        cycle store_basis
     end if
  end do store_basis

  !close the file of the basis definition
  close(35)

!!$  !print the found values
!!$  do ityp=1,ntypes
!!$     do ishell=1,nshell(ityp)
!!$        print *,'ityp=',ityp,'ishell=',ishell,'l=',nam(ishell,ityp),'ndoc=',ndoc(ishell,ityp)
!!$        do ipg=1,ndoc(ishell,ityp)
!!$           print *,'expo=',expo(ipg,ishell,ityp),'coeff=',contcoeff(ipg,ishell,ityp)
!!$        end do
!!$     end do
!!$  end do

  mmax=2*lmax+1
  !now read the coefficients of the gaussian converged orbitals
  open(unit=36,file='gaucoeff.dat',action='read')
  !here there is the orbital label, for the moment it is assumed to vary between 1 and 4
  allocate(ctmp(10),stat=i_stat)
  call memocc(i_stat,product(shape(ctmp))*kind(ctmp),'ctmp','gautowav')
  allocate(iorbtmp(10),stat=i_stat)
  call memocc(i_stat,product(shape(iorbtmp))*kind(iorbtmp),'iorbtmp','gautowav')
  allocate(cimu(mmax,nbx,nat,norb),stat=i_stat)
  call memocc(i_stat,product(shape(cimu))*kind(cimu),'cimu','gautowav')

  read(36,*)
  read_line1: do
     read(36,'(a100)')line
     !analyse how many orbitals are contained in a given line
     read_orbitals1: do ipar=10,1,-1
        read(line,*,iostat=i_stat)(iorbtmp(i),i=1,ipar)
        if (i_stat==0) then
           read(line,*)nst
           exit read_line1
        end if
     end do read_orbitals1
  end do read_line1
  nend=nst+ipar-1

!!$  nst=1
!!$  nend=4
  jat=0
  ishell=1
  jbas=0
  iat=nat
  !now read the data for assignating the coefficients
  store_coeff: do
     read(36,'(a100)')line
     !choose between different cases
     read(line,*,iostat=i_stat)ibas,iat,symbol,string,(ctmp(iorb),iorb=1,nend-nst+1)
     if (i_stat==0) then
        !print *,line,nst,nend
        if (jat==iat) then
           jbas=jbas+1
           if (jbas > 2*nam(ishell,iatype(iat))+1) then
              jbas=1
              ishell=ishell+1
              if (ishell > nshell(iatype(iat))) then
                 if (iproc==0) write(*,'(1x,a,i0,a)')&
                      'Problem in the gaucoeff.dat file, the number of shells of atom ',iat ,&
                      ' is incoherent'
                 stop
              end if
           end if
        else
           jbas=1
           ishell=1
        end if
       symbol=trim(string)
       do iorb=nst,nend
          cimu(jbas+myshift(symbol),ishell,iat,iorb)=ctmp(iorb-nst+1)
       end do
       jat=iat
       if (jbas==2*nam(ishell,iatype(iat))+1 .and. ishell==nshell(iatype(iat))&
            .and. iat==nat .and. nend==norb) then
          exit store_coeff
       else
          cycle store_coeff
       end if
     end if

     read_orbitals: do ipar=10,1,-1
        read(line,*,iostat=i_stat)(iorbtmp(i),i=1,ipar)
        if (i_stat==0) then
           read(line,*)nst
           nend=nst+ipar-1
           if (jat/=nat) then
              if (iproc==0)write(*,'(1x,a,i0,a)')&
                   'Problem in the gaucoeff.dat file, only ',iat ,' atoms processed'
              stop
           else
              cycle store_coeff
           end if
        end if
     end do read_orbitals

  end do store_coeff
  close(36)

!!$  !print the found values
!!$  do iat=1,nat
!!$     ityp=iatype(iat)
!!$     do ishell=1,nshell(ityp)
!!$        do jbas=1,2*nam(ishell,ityp)+1
!!$           print *,iat,ishell,nam(ishell,ityp),jbas,(cimu(jbas,ishell,iat,iorb),iorb=1,norb)
!!$        end do
!!$     end do
!!$  end do

  i_all=-product(shape(ctmp))*kind(ctmp)
  deallocate(ctmp,stat=i_stat)
  call memocc(i_stat,i_all,'ctmp','gautowav')
  i_all=-product(shape(iorbtmp))*kind(iorbtmp)
  deallocate(iorbtmp,stat=i_stat)
  call memocc(i_stat,i_all,'iorbtmp','gautowav')


  !now apply this basis set information to construct the wavelets wavefunctions
  allocate(lx(nterm_max),stat=i_stat)
  call memocc(i_stat,product(shape(lx))*kind(lx),'lx','gautowav')
  allocate(ly(nterm_max),stat=i_stat)
  call memocc(i_stat,product(shape(ly))*kind(ly),'ly','gautowav')
  allocate(lz(nterm_max),stat=i_stat)
  call memocc(i_stat,product(shape(lz))*kind(lz),'lz','gautowav')
  allocate(fac_arr(nterm_max),stat=i_stat)
  call memocc(i_stat,product(shape(fac_arr))*kind(fac_arr),'fac_arr','gautowav')

  allocate(psiatn(ngx),stat=i_stat)
  call memocc(i_stat,product(shape(psiatn))*kind(psiatn),'psiatn','gautowav')
  allocate(xp(ngx),stat=i_stat)
  call memocc(i_stat,product(shape(xp))*kind(xp),'xp','gautowav')
  allocate(tpsi(nvctr_c+7*nvctr_f),stat=i_stat)
  call memocc(i_stat,product(shape(tpsi))*kind(tpsi),'tpsi','gautowav')


  !initialize the wavefunction
  psi(:,:)=0.d0

  

  !this can be changed to be passed only once to all the gaussian basis
  eks=0.d0
  !loop over the atoms
  do iat=1,nat
     ityp=iatype(iat)
     rx=rxyz(1,iat)
     ry=rxyz(2,iat)
     rz=rxyz(3,iat)
     !loop over the number of shells of the atom type
     do ishell=1,nshell(ityp)
        !the degree of contraction of the basis function
        !is the same as the ng value of the createAtomicOrbitals routine
        ng=ndoc(ishell,ityp)
        !angular momentum of the basis set(shifted for compatibility with BigDFT routines
        l=nam(ishell,ityp)+1
        !amplitude coefficients (contraction coefficients of the basis times
        !the amplitude of this basis in that orbital)
        !exponents for the gaussian expansion adapted following the convention
        !of the routine gauss_to_daub
        do ig=1,ng
           psiatn(ig)=contcoeff(ig,ishell,ityp)
           xp(ig)=sqrt(0.5d0/expo(ig,ishell,ityp))
        end do
        call atomkin(l-1,ng,xp,psiatn,psiatn,ek)
        !multiply the values of the gaussian contraction times the orbital coefficient
        do m=1,2*l-1
           call calc_coeff_inguess(l,m,nterm_max,nterm,lx,ly,lz,fac_arr)
!!$           !this kinetic energy is not reliable
!!$           eks=eks+ek*occup(iorb)*cimu(m,ishell,iat,iorb)
           call crtonewave(n1,n2,n3,ng,nterm,lx,ly,lz,fac_arr,xp,psiatn,&
                rx,ry,rz,hgrid,0,n1,0,n2,0,n3,nfl1,nfu1,nfl2,nfu2,nfl3,nfu3,  & 
                nseg_c,nvctr_c,keyg,keyv,nseg_f,nvctr_f,&
                keyg(1,nseg_c+1),keyv(nseg_c+1),&
                tpsi(1),tpsi(nvctr_c+1))
           !sum the result inside the orbital wavefunction
           !loop over the orbitals
           do iorb=1,norb
              if (myorbital(iorb,norb,iproc,nproc)) then
                 jorb=iorb-iproc*norbp
                 do i=1,nvctr_c+7*nvctr_f
                    psi(i,jorb)=psi(i,jorb)+cimu(m,ishell,iat,iorb)*tpsi(i)
                 end do
              end if
           end do
        end do
     end do
     if (iproc == 0) then
        write(*,'(a)',advance='no') &
             repeat('.',(iat*40)/nat-((iat-1)*40)/nat)
     end if
  end do
  !renormalize the orbitals
  !write the norm of the orbital
  do iorb=1,norb
     if (myorbital(iorb,norb,iproc,nproc)) then
        jorb=iorb-iproc*norbp
        call wnrm(nvctr_c,nvctr_f,psi(1,jorb),psi(nvctr_c+1,jorb),scpr) 
        call wscal(nvctr_c,nvctr_f,1.d0/sqrt(scpr),psi(1,jorb),psi(nvctr_c+1,jorb))
        print *,'norm of orbital ',iorb,scpr
     end if
  end do

  !now we have to evaluate the eigenvalues of this hamiltonian

  i_all=-product(shape(tpsi))*kind(tpsi)
  deallocate(tpsi,stat=i_stat)
  call memocc(i_stat,i_all,'tpsi','gautowav')

  i_all=-product(shape(nshell))*kind(nshell)
  deallocate(nshell,stat=i_stat)
  call memocc(i_stat,i_all,'nshell','gautowav')
  i_all=-product(shape(nam))*kind(nam)
  deallocate(nam,stat=i_stat)
  call memocc(i_stat,i_all,'nam','gautowav')
  i_all=-product(shape(ndoc))*kind(ndoc)
  deallocate(ndoc,stat=i_stat)
  call memocc(i_stat,i_all,'ndoc','gautowav')
  i_all=-product(shape(contcoeff))*kind(contcoeff)
  deallocate(contcoeff,stat=i_stat)
  call memocc(i_stat,i_all,'contcoeff','gautowav')
  i_all=-product(shape(expo))*kind(expo)
  deallocate(expo,stat=i_stat)
  call memocc(i_stat,i_all,'expo','gautowav')
  i_all=-product(shape(cimu))*kind(cimu)
  deallocate(cimu,stat=i_stat)
  call memocc(i_stat,i_all,'cimu','gautowav')

  i_all=-product(shape(lx))*kind(lx)
  deallocate(lx,stat=i_stat)
  call memocc(i_stat,i_all,'lx','gautowav')
  i_all=-product(shape(ly))*kind(ly)
  deallocate(ly,stat=i_stat)
  call memocc(i_stat,i_all,'ly','gautowav')
  i_all=-product(shape(lz))*kind(lz)
  deallocate(lz,stat=i_stat)
  call memocc(i_stat,i_all,'lz','gautowav')
  i_all=-product(shape(fac_arr))*kind(fac_arr)
  deallocate(fac_arr,stat=i_stat)
  call memocc(i_stat,i_all,'fac_arr','gautowav')

  i_all=-product(shape(xp))*kind(xp)
  deallocate(xp,stat=i_stat)
  call memocc(i_stat,i_all,'xp','gautowav')
  i_all=-product(shape(psiatn))*kind(psiatn)
  deallocate(psiatn,stat=i_stat)
  call memocc(i_stat,i_all,'psiatn','gautowav')


end subroutine gautowav

subroutine atomkin(l,ng,xp,psiat,psiatn,ek)
  ! calculates the kinetic energy of an atomic wavefunction expressed in Gaussians
  ! the output psiatn is a normalized version of psiat
  implicit real*8 (a-h,o-z)
  dimension xp(ng),psiat(ng),psiatn(ng)

  !        gml=.5d0*gamma(.5d0+l)
  gml = 0.d0
  if (l.eq.0) then 
     gml=0.88622692545275801365d0
  else if (l.eq.1) then 
     gml=0.44311346272637900682d0
  else if (l.eq.2) then 
     gml=0.66467019408956851024d0
  else if (l.eq.3) then 
     gml=1.6616754852239212756d0
  else
     stop 'atomkin'
  endif

  ek=0.d0
  tt=0.d0
  do i=1,ng
     xpi=.5d0/xp(i)**2
     do j=1,ng
        xpj=.5d0/xp(j)**2
        d=xpi+xpj
        sxp=1.d0/d
        const=gml*sqrt(sxp)**(2*l+1)
        ! kinetic energy  matrix element hij
        hij=.5d0*const*sxp**2* ( 3.d0*xpi*xpj +                  &
             l*(6.d0*xpi*xpj-xpi**2-xpj**2) -        &
             l**2*(xpi-xpj)**2  ) + .5d0*l*(l+1.d0)*const
        sij=const*sxp*(l+.5d0)
        ek=ek+hij*psiat(i)*psiat(j)
        tt=tt+sij*psiat(i)*psiat(j)
     enddo
  enddo

  !commented out, to be reinserted
  !if (abs(tt-1.d0).gt.1.d-2) write(*,*) 'presumably wrong inguess data',l,tt
  ! energy expectation value
  ek=ek/tt
  !write(*,*) 'ek=',ek,tt,l,ng
  ! scale atomic wavefunction
  tt=sqrt(1.d0/tt)
!!$        if (l.eq.0) then  ! multiply with 1/sqrt(4*pi)
!!$        tt=tt*0.28209479177387814347d0
!!$        else if (l.eq.1) then  ! multiply with sqrt(3/(4*pi))
!!$        tt=tt*0.48860251190291992159d0
!!$        !decide the value of the normalization to be used
!!$        endif
  do i=1,ng
     psiatn(i)=psiat(i)*tt
  enddo

  return
end subroutine atomkin

!calculate the shift between the spherical harmonics of CP2K and the one of BigDFT
function myshift(symbol)
  implicit none
  character(len=5), intent(in) :: symbol
  integer :: myshift
  if (symbol(2:2)=='s') then
     myshift=0
  else if (symbol(2:2)=='p') then
     if( symbol(3:3)=='y') then
        myshift=1
     else if( symbol(3:3)=='z') then
        myshift=1
     else if( symbol(3:3)=='x') then
        myshift=-2
     end if
  else if ( symbol(2:2)=='d') then
     if( symbol(3:4)=='-2') then
        myshift=2
     else if( symbol(3:4)=='-1') then
        myshift=-1
     else if( symbol(3:3)=='0') then
        myshift=2
     else if( symbol(3:4)=='+1') then
        myshift=-2
     else if( symbol(3:4)=='+2') then
        myshift=-1
     end if
  else if ( symbol(2:2)=='f') then
     if( symbol(3:4)=='-3') then
        myshift=4
     else if( symbol(3:4)=='-2') then
        myshift=5
     else if( symbol(3:4)=='-1') then
        myshift=-1
     else if( symbol(3:3)=='0') then
        myshift=-1
     else if( symbol(3:4)=='+1') then
        myshift=-4
     else if( symbol(3:4)=='+2') then
        myshift=0
     else if( symbol(3:4)=='+3') then
        myshift=-3
     end if
  else if ( symbol(2:2)=='g') then
     write(*,'(1x,a)')'the orbitals of type g are not yet implemented in BigDFT'
     stop
  end if

end function myshift

logical function myorbital(iorb,norbe,iproc,nproc)
  implicit real*8 (a-h,o-z)
  parameter(eps_mach=1.d-12)

  tt=dble(norbe)/dble(nproc)
  norbep=int((1.d0-eps_mach*tt) + tt)
  if (iorb .ge. iproc*norbep+1 .and. iorb .le. min((iproc+1)*norbep,norbe)) then
     myorbital=.true.
  else
     myorbital=.false.
  endif

  return
end function myorbital

subroutine crtonewave(n1,n2,n3,nterm,ntp,lx,ly,lz,fac_arr,xp,psiat,rx,ry,rz,hgrid, & 
     nl1_c,nu1_c,nl2_c,nu2_c,nl3_c,nu3_c,nl1_f,nu1_f,nl2_f,nu2_f,nl3_f,nu3_f,  & 
     nseg_c,mvctr_c,keyg_c,keyv_c,nseg_f,mvctr_f,keyg_f,keyv_f,psi_c,psi_f)
  ! returns an input guess orbital that is a Gaussian centered at a Wannier center
  ! exp (-1/(2*gau_a^2) *((x-cntrx)^2 + (y-cntry)^2 + (z-cntrz)^2 ))
  ! in the arrays psi_c, psi_f
  implicit real*8 (a-h,o-z)
  parameter(nw=16000)
  dimension xp(nterm),psiat(nterm),fac_arr(ntp)
  dimension lx(ntp),ly(ntp),lz(ntp)
  dimension keyg_c(2,nseg_c),keyv_c(nseg_c),keyg_f(2,nseg_f),keyv_f(nseg_f)
  dimension psi_c(mvctr_c),psi_f(7,mvctr_f)
  real*8, allocatable, dimension(:,:) :: wprojx, wprojy, wprojz
  real*8, allocatable, dimension(:,:) :: work
  real*8, allocatable :: psig_c(:,:,:), psig_f(:,:,:,:)

  allocate(wprojx(0:n1,2),stat=i_stat)
  call memocc(i_stat,product(shape(wprojx))*kind(wprojx),'wprojx','crtonewave')
  allocate(wprojy(0:n2,2),stat=i_stat)
  call memocc(i_stat,product(shape(wprojy))*kind(wprojy),'wprojy','crtonewave')
  allocate(wprojz(0:n3,2),stat=i_stat)
  call memocc(i_stat,product(shape(wprojz))*kind(wprojz),'wprojz','crtonewave')
  allocate(work(0:nw,2),stat=i_stat)
  call memocc(i_stat,product(shape(work))*kind(work),'work','crtonewave')
  allocate(psig_c(nl1_c:nu1_c,nl2_c:nu2_c,nl3_c:nu3_c),stat=i_stat)
  call memocc(i_stat,product(shape(psig_c))*kind(psig_c),'psig_c','crtonewave')
  allocate(psig_f(7,nl1_f:nu1_f,nl2_f:nu2_f,nl3_f:nu3_f),stat=i_stat)
  call memocc(i_stat,product(shape(psig_f))*kind(psig_f),'psig_f','crtonewave')

  iterm=1
  itp=1
  gau_a=xp(iterm)
  n_gau=lx(itp)
  CALL GAUSS_TO_DAUB(hgrid,fac_arr(itp),rx,gau_a,n_gau,n1,ml1,mu1,wprojx(0,1),te,work,nw)
  n_gau=ly(itp)
  CALL GAUSS_TO_DAUB(hgrid,1.d0,ry,gau_a,n_gau,n2,ml2,mu2,wprojy(0,1),te,work,nw)
  n_gau=lz(itp)
  CALL GAUSS_TO_DAUB(hgrid,psiat(iterm),rz,gau_a,n_gau,n3,ml3,mu3,wprojz(0,1),te,work,nw)

  ! First term: coarse projector components
  do i3=nl3_c,nu3_c
     do i2=nl2_c,nu2_c
        do i1=nl1_c,nu1_c
           psig_c(i1,i2,i3)=wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,1)
        enddo
     enddo
  enddo

  ! First term: fine projector components
  do i3=nl3_f,nu3_f
     do i2=nl2_f,nu2_f
        do i1=nl1_f,nu1_f
           psig_f(1,i1,i2,i3)=wprojx(i1,2)*wprojy(i2,1)*wprojz(i3,1)
           psig_f(2,i1,i2,i3)=wprojx(i1,1)*wprojy(i2,2)*wprojz(i3,1)
           psig_f(3,i1,i2,i3)=wprojx(i1,2)*wprojy(i2,2)*wprojz(i3,1)
           psig_f(4,i1,i2,i3)=wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,2)
           psig_f(5,i1,i2,i3)=wprojx(i1,2)*wprojy(i2,1)*wprojz(i3,2)
           psig_f(6,i1,i2,i3)=wprojx(i1,1)*wprojy(i2,2)*wprojz(i3,2)
           psig_f(7,i1,i2,i3)=wprojx(i1,2)*wprojy(i2,2)*wprojz(i3,2)
        enddo
     enddo
  enddo

  do iterm=2,nterm
     gau_a=xp(iterm)
     n_gau=lx(itp)
     CALL GAUSS_TO_DAUB(hgrid,fac_arr(itp),rx,gau_a,n_gau,n1,ml1,mu1,wprojx(0,1),te,work,nw)
     n_gau=ly(itp)
     CALL GAUSS_TO_DAUB(hgrid,1.d0,ry,gau_a,n_gau,n2,ml2,mu2,wprojy(0,1),te,work,nw)
     n_gau=lz(itp)
     CALL GAUSS_TO_DAUB(hgrid,psiat(iterm),rz,gau_a,n_gau,n3,ml3,mu3,wprojz(0,1),te,work,nw)

     ! First term: coarse projector components
     do i3=nl3_c,nu3_c
        do i2=nl2_c,nu2_c
           do i1=nl1_c,nu1_c
              psig_c(i1,i2,i3)=psig_c(i1,i2,i3)+wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,1)
           enddo
        enddo
     enddo

     ! First term: fine projector components
     do i3=nl3_f,nu3_f
        do i2=nl2_f,nu2_f
           do i1=nl1_f,nu1_f
              psig_f(1,i1,i2,i3)=psig_f(1,i1,i2,i3)+wprojx(i1,2)*wprojy(i2,1)*wprojz(i3,1)
              psig_f(2,i1,i2,i3)=psig_f(2,i1,i2,i3)+wprojx(i1,1)*wprojy(i2,2)*wprojz(i3,1)
              psig_f(3,i1,i2,i3)=psig_f(3,i1,i2,i3)+wprojx(i1,2)*wprojy(i2,2)*wprojz(i3,1)
              psig_f(4,i1,i2,i3)=psig_f(4,i1,i2,i3)+wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,2)
              psig_f(5,i1,i2,i3)=psig_f(5,i1,i2,i3)+wprojx(i1,2)*wprojy(i2,1)*wprojz(i3,2)
              psig_f(6,i1,i2,i3)=psig_f(6,i1,i2,i3)+wprojx(i1,1)*wprojy(i2,2)*wprojz(i3,2)
              psig_f(7,i1,i2,i3)=psig_f(7,i1,i2,i3)+wprojx(i1,2)*wprojy(i2,2)*wprojz(i3,2)
           enddo
        enddo
     enddo

  end do

  do itp=2,ntp

     do iterm=1,nterm
        gau_a=xp(iterm)
        n_gau=lx(itp)
        CALL GAUSS_TO_DAUB(hgrid,fac_arr(itp),rx,gau_a,n_gau,n1,ml1,mu1,wprojx(0,1),te,work,nw)
        n_gau=ly(itp)
        CALL GAUSS_TO_DAUB(hgrid,1.d0,ry,gau_a,n_gau,n2,ml2,mu2,wprojy(0,1),te,work,nw)
        n_gau=lz(itp)
        CALL GAUSS_TO_DAUB(hgrid,psiat(iterm),rz,gau_a,n_gau,n3,ml3,mu3,wprojz(0,1),te,work,nw)

        ! First term: coarse projector components
        do i3=nl3_c,nu3_c
           do i2=nl2_c,nu2_c
              do i1=nl1_c,nu1_c
                 psig_c(i1,i2,i3)=psig_c(i1,i2,i3)+wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,1)
              enddo
           enddo
        enddo

        ! First term: fine projector components
        do i3=nl3_f,nu3_f
           do i2=nl2_f,nu2_f
              do i1=nl1_f,nu1_f
                 psig_f(1,i1,i2,i3)=psig_f(1,i1,i2,i3)+wprojx(i1,2)*wprojy(i2,1)*wprojz(i3,1)
                 psig_f(2,i1,i2,i3)=psig_f(2,i1,i2,i3)+wprojx(i1,1)*wprojy(i2,2)*wprojz(i3,1)
                 psig_f(3,i1,i2,i3)=psig_f(3,i1,i2,i3)+wprojx(i1,2)*wprojy(i2,2)*wprojz(i3,1)
                 psig_f(4,i1,i2,i3)=psig_f(4,i1,i2,i3)+wprojx(i1,1)*wprojy(i2,1)*wprojz(i3,2)
                 psig_f(5,i1,i2,i3)=psig_f(5,i1,i2,i3)+wprojx(i1,2)*wprojy(i2,1)*wprojz(i3,2)
                 psig_f(6,i1,i2,i3)=psig_f(6,i1,i2,i3)+wprojx(i1,1)*wprojy(i2,2)*wprojz(i3,2)
                 psig_f(7,i1,i2,i3)=psig_f(7,i1,i2,i3)+wprojx(i1,2)*wprojy(i2,2)*wprojz(i3,2)
              enddo
           enddo
        enddo

     end do


  end do


  !wavefunction compression
  ! coarse part
  do iseg=1,nseg_c
     jj=keyv_c(iseg)
     j0=keyg_c(1,iseg)
     j1=keyg_c(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psi_c(i-i0+jj)=psig_c(i,i2,i3)
     enddo
  enddo

  ! fine part
  do iseg=1,nseg_f
     jj=keyv_f(iseg)
     j0=keyg_f(1,iseg)
     j1=keyg_f(2,iseg)
     ii=j0-1
     i3=ii/((n1+1)*(n2+1))
     ii=ii-i3*(n1+1)*(n2+1)
     i2=ii/(n1+1)
     i0=ii-i2*(n1+1)
     i1=i0+j1-j0
     do i=i0,i1
        psi_f(1,i-i0+jj)=psig_f(1,i,i2,i3)
        psi_f(2,i-i0+jj)=psig_f(2,i,i2,i3)
        psi_f(3,i-i0+jj)=psig_f(3,i,i2,i3)
        psi_f(4,i-i0+jj)=psig_f(4,i,i2,i3)
        psi_f(5,i-i0+jj)=psig_f(5,i,i2,i3)
        psi_f(6,i-i0+jj)=psig_f(6,i,i2,i3)
        psi_f(7,i-i0+jj)=psig_f(7,i,i2,i3)
     enddo
  enddo

  i_all=-product(shape(wprojx))*kind(wprojx)
  deallocate(wprojx,stat=i_stat)
  call memocc(i_stat,i_all,'wprojx','crtonewave')
  i_all=-product(shape(wprojy))*kind(wprojy)
  deallocate(wprojy,stat=i_stat)
  call memocc(i_stat,i_all,'wprojy','crtonewave')
  i_all=-product(shape(wprojz))*kind(wprojz)
  deallocate(wprojz,stat=i_stat)
  call memocc(i_stat,i_all,'wprojz','crtonewave')
  i_all=-product(shape(work))*kind(work)
  deallocate(work,stat=i_stat)
  call memocc(i_stat,i_all,'work','crtonewave')
  i_all=-product(shape(psig_c))*kind(psig_c)
  deallocate(psig_c,stat=i_stat)
  call memocc(i_stat,i_all,'psig_c','crtonewave')
  i_all=-product(shape(psig_f))*kind(psig_f)
  deallocate(psig_f,stat=i_stat)
  call memocc(i_stat,i_all,'psig_f','crtonewave')

  return
END SUBROUTINE crtonewave


subroutine calc_coeff_inguess(l,m,nterm_max,nterm,lx,ly,lz,fac_arr)

  implicit none
  integer, intent(in) :: l,m,nterm_max
  integer, intent(out) :: nterm
  integer, dimension(nterm_max), intent(out) :: lx,ly,lz
  real(kind=8), dimension(nterm_max), intent(out) :: fac_arr

  if (l.eq.1 .and. m.eq.1) then
     nterm=1
     lx(1)=0 ; ly(1)=0 ; lz(1)=0
     fac_arr(1)=0.28209479177387814347d0

  else if (l.eq.2  .and. m.eq.1) then
     nterm=1
     lx(1)=1 ; ly(1)=0 ; lz(1)=0
     fac_arr(1)=0.48860251190291992159d0
  else if (l.eq.2  .and. m.eq.2) then
     nterm=1
     lx(1)=0 ; ly(1)=1 ; lz(1)=0
     fac_arr(1)=0.48860251190291992159d0
  else if (l.eq.2  .and. m.eq.3) then
     nterm=1
     lx(1)=0 ; ly(1)=0 ; lz(1)=1
     fac_arr(1)=0.48860251190291992159d0

  else if (l.eq.3  .and. m.eq.1) then
     nterm=1
     lx(1)=0 ; ly(1)=1 ; lz(1)=1
     fac_arr(1)=1.092548430592079d0
  else if (l.eq.3  .and. m.eq.2) then
     nterm=1
     lx(1)=1 ; ly(1)=0 ; lz(1)=1
     fac_arr(1)=1.092548430592079d0
  else if (l.eq.3  .and. m.eq.3) then
     nterm=1
     lx(1)=1 ; ly(1)=1 ; lz(1)=0
     fac_arr(1)=1.092548430592079d0
  else if (l.eq.3  .and. m.eq.4) then
     nterm=2
     lx(1)=2 ; ly(1)=0 ; lz(1)=0
     lx(2)=0 ; ly(2)=2 ; lz(2)=0
     fac_arr(1)=0.5462742152960396d0
     fac_arr(2)=-0.5462742152960396d0
  else if (l.eq.3  .and. m.eq.5) then 
     nterm=3
     lx(1)=2 ; ly(1)=0 ; lz(1)=0
     lx(2)=0 ; ly(2)=2 ; lz(2)=0
     lx(3)=0 ; ly(3)=0 ; lz(3)=2
     fac_arr(1)=-0.3153915652525201d0
     fac_arr(2)=-0.3153915652525201d0
     fac_arr(3)=2.d0*0.3153915652525201d0

  else if (l.eq.4  .and. m.eq.1) then
     nterm=3
     lx(1)=3 ; ly(1)=0 ; lz(1)=0
     lx(2)=1 ; ly(2)=2 ; lz(2)=0
     lx(3)=1 ; ly(3)=0 ; lz(3)=2
     fac_arr(1)=0.4570457994644658d0
     fac_arr(2)=0.4570457994644658d0
     fac_arr(3)=-4.d0*0.4570457994644658d0
  else if (l.eq.4  .and. m.eq.2) then
     nterm=3
     lx(1)=2 ; ly(1)=1 ; lz(1)=0
     lx(2)=0 ; ly(2)=3 ; lz(2)=0
     lx(3)=0 ; ly(3)=1 ; lz(3)=2
     fac_arr(1)=0.4570457994644658d0
     fac_arr(2)=0.4570457994644658d0
     fac_arr(3)=-4.d0*0.4570457994644658d0
  else if (l.eq.4  .and. m.eq.3) then
     nterm=3
     lx(1)=2 ; ly(1)=0 ; lz(1)=1
     lx(2)=0 ; ly(2)=2 ; lz(2)=1
     lx(3)=0 ; ly(3)=0 ; lz(3)=3
     fac_arr(1)=3.d0*0.3731763325901154d0
     fac_arr(2)=3.d0*0.3731763325901154d0
     fac_arr(3)=-2.d0*0.3731763325901154d0
  else if (l.eq.4  .and. m.eq.4) then
     nterm=2
     lx(1)=3 ; ly(1)=0 ; lz(1)=0
     lx(2)=1 ; ly(2)=2 ; lz(2)=0
     fac_arr(1)=0.5900435899266436d0
     fac_arr(2)=-3.d0*0.5900435899266436d0
  else if (l.eq.4  .and. m.eq.5) then
     nterm=2
     lx(1)=2 ; ly(1)=1 ; lz(1)=0
     lx(2)=0 ; ly(2)=3 ; lz(2)=0
     fac_arr(1)=-3.d0*0.5900435899266436d0
     fac_arr(2)=0.5900435899266436d0
  else if (l.eq.4  .and. m.eq.6) then
     nterm=2
     lx(1)=2 ; ly(1)=0 ; lz(1)=1
     lx(2)=0 ; ly(2)=2 ; lz(2)=1
     fac_arr(1)=1.445305721320277d0
     fac_arr(2)=-1.445305721320277d0
  else if (l.eq.4  .and. m.eq.7) then
     nterm=1
     lx(1)=1 ; ly(1)=1 ; lz(1)=1
     fac_arr(1)=2.890611442640554d0
  else
     stop 'input guess format error'
  endif

END SUBROUTINE calc_coeff_inguess

