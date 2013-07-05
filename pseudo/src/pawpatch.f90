!> @file
!! Paw generation (?) in pseudo program
!! @author
!!    Alex Willand, under the supervision of Stefan Goedecker
!!    gpu accelerated routines by Raffael Widmer
!!    parts of this program were based on the fitting program by matthias krack
!!    http://cvs.berlios.de/cgi-bin/viewcvs.cgi/cp2k/potentials/goedecker/pseudo/v2.2/
!!
!!    Copyright (C) 2010-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


subroutine pawpatch(energ,verbose,maxdim,pp,penal,&
     noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,pol,nsmx,&
     no,lo,so,ev,crcov,dcrcov,ddcrcov,norb,&
     occup,aeval,chrg,dhrg,ehrg,res,wght,&
     wfnode,psir0,wghtp0,&
     rcov,rprb,rcore,zcore,znuc,zion,rloc,gpot,r_l,hsep,&
     vh,xp,rmt,rmtg,ud,nint,ng_fit,ngmx,psi,&
     avgl1,avgl2,avgl3,ortprj,litprj,igrad,rae, &
     iproc,nproc,wghtconf,wghtexci,wghtsoft,wghtrad,wghthij,&
     nhgrid,hgridmin,hgridmax, nhpow,ampl,crmult,frmult,&
     excitAE,ntime,iter,itertot,penref,time,ngrid_fit,&
     nconfpaw, npawl, nchannelspaw, ispp, pawstatom , pawstN, pawstL, pawstP ,&
     pawrcovfact)
  
  implicit none
  !! implicit real*8 (a-h,o-z)
  integer, parameter :: gp=kind(1.0d0) 
  logical avgl1,avgl2,avgl3,ortprj,litprj,igrad
  real(8) pp(maxdim),so(norb),ev(norb),crcov(norb),&
       dcrcov(norb),ddcrcov(norb),occup(noccmx,lmx,nsmx),aeval(noccmx,lmx,nsmx),&
       chrg(noccmx,lmx,nsmx),dhrg(noccmx,lmx,nsmx),&
       ehrg(noccmx,lmx,nsmx),res(noccmx,lmx,nsmx),&
       wght(noccmx,lmx,nsmx,8),&
       wfnode(noccmx,lmx,nsmx,3),&
       gpot(*),r_l(*),hsep(6,lpmx,nsmx),&
       vh(*),xp(*),rmt(*),rmtg(*),ud(*),psi(*),&
       rae(*),&
       time(3), penal, psir0,wghtp0,rcov,&
       rprb,rcore,zcore,znuc,zion,rloc,&
       wghtexci,wghtsoft,wghtrad,wghthij,&
       hgridmin,hgridmax, ampl,crmult,frmult,&
       excitAE,iter,itertot,penref, wghtconf
  
  integer no(norb),lo(norb),nconfpaw, npawl, nchannelspaw , maxdim,&
       noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx,nint,ng_fit,ngmx,iproc,&
       nproc,nhgrid,nhpow,ntime, norb, ngrid_fit, j,   pawstN, pawstL, pawstP
  
  real(8) pawrcovfact

  character(len=1) :: ispp
  logical :: energ, verbose, pol
  character(len=30) :: plotfile
  real(8), pointer :: atom_potential_fit(:)
  real(8), pointer :: statom_potential(:)
  real(8) :: rdum
  
  integer :: Npaw, ng
  integer Ngrid, Ngrid_box, Ngrid_biggerbox, iovrsmpl,  Ngrid_box_larger
  real(8) boxradius, biggerboxradius, a,b, EMAX
  real(8), pointer :: rgrid(:), yp(:), ypp(:), w(:), aepot(:), aepot_p(:), aepot_pp(:), &
       rgrid_ab(:), aepot_cent(:), staepot(:), rw(:),rd(:)
  real(8) a1,b1,an,bn
  integer ierr, isx
  integer LPaw, n, Nsol
  integer igrid
  real(8), pointer :: psi_initial_copy(:),psi_initial(:), dumpsi_p(:)
  real(8) dum_energy
  character(1000) filename
  character(len=125) :: pawstatom
  real(gp) , pointer ::  psigrid(:,:), Egrid(:), nonloc(:)
  real(gp) , pointer ::  psigrid_pseudo(:,:), Egrid_pseudo(:)
  real(gp) , pointer ::  psipsigrid_pseudo(:,:)
  real(gp) , pointer ::  psigrid_bigger(:,:)
  
  real(gp) , pointer ::  expo(:)
  real(gp), pointer::PAWpatch_matrix(:,:)
  real(gp) fourpi
  logical dump_functions
  integer real_start, l, iocc,i 
  integer iout, outunits(2), uout
  include 'mpif.h'
  
  dump_functions= .true.

  
  !! for verbose output 
  open(unit=38, file="pawpatch.verbose" )
  outunits(1)=6
  outunits(2)=38

  do iout=1,2
     uout=outunits(iout)
     write(uout,*) "            "
     write(uout,*) "---PAWPATCH  parameters--- "
     write(uout,*) " nconfpaw     ", nconfpaw
     write(uout,*) " npawl        ", npawl
     write(uout,*) " nchannelspaw ", nchannelspaw
     write(uout,*) " pawstatom    ", trim(pawstatom)
     write(uout,*) " pawstN       ", pawstN
     write(uout,*) " pawstL       ", pawstL 
     write(uout,*) " pawstP       ", pawstP
     write(uout,*) " pawrcovfact  ", pawrcovfact
    write(uout,*) "            "
  end do


  if (nspin/=1) then
     ! stop "pawpatch can be built only in the nspin=1 case   " 
     write(6,* )  , "WARNING : pawpatch works in the nspin=1 case, now moving to this case modifying occup and hsep  "      
     write(38,*)    "WARNING : pawpatch works in the nspin=1 case, now moving to this case modifying occup and hsep  "      
     !! reducing occup for nspin=1
     do l=0,lmax
        do iocc=1,noccmax
           occup(iocc,l+1,1) =occup(iocc,l+1,1)+ occup(iocc,l+1,2)
        end do
     end do
     
     do i=1,6
        do l=1, lmax
           hsep(i,l,1)= ( hsep(i,l,1)*l +   hsep(i,l,2)*(l-1)  ) /(2*l-1)
        end do
     end do
  endif


  if(iproc/=0) return
  allocate(atom_potential_fit(ngrid_fit))
  write(plotfile, '(a,i0,a)') 'ae.pot.conf.',nconfpaw ,'.plt'
  open(unit=37,file=trim(plotfile),status='unknown')
  do j=1,ngrid_fit
     read(37, *)   rdum, atom_potential_fit(j)
     if( abs(rdum-rae(j))>1.0e-6) then
        STOP "rgrid not corresponding "
     end if
  end do
  close(37)
  
  allocate(statom_potential(ngrid_fit))
  if ( trim(pawstatom)/='') then
     open(unit=37,file=trim(pawstatom),status='unknown')
     do j=1,ngrid_fit
        read(37, *)   rdum, statom_potential(j)
        if( abs(rdum-rae(j))>1.0e-6) then
           STOP "rgrid in statom  not corresponding  with the one in ae.pot.conf "
        end if
     end do
     close(37)
  endif

  !! partly correct for r(1) to be 0 and potential(r=0) meaningless
  atom_potential_fit(1)=2*atom_potential_fit(2)
  statom_potential  (1)=2*statom_potential  (2)

  ng  = 30
  !! noccmax = 5 
  !! print *, "  NOCCMAX as set by pseudo is  ", noccmax
  !! lmax=3
  !! print *,"LMAX as set by pseudo is ", lmax 
  
  Nsol=100
  
  Npaw= nchannelspaw
  
  Ngrid=10000                      !! roughly the number of point of the 
  !! oversampled grid
  
  boxradius=rcov *pawrcovfact                  !! this should be found in the original grid
  biggerboxradius = 1.5_8 * rcov   !! this is an approximative goal
  
  Ngrid_box=1
  Ngrid_box_larger=1
  Ngrid_biggerbox=1

  
  do j=1,ngrid_fit
     if( abs(rae(j) -boxradius) <  abs(rae(Ngrid_box) -boxradius)) Ngrid_box = j
     ! if( abs(rae(j) -boxradius*2) <  abs(rae(Ngrid_box_larger) -boxradius*2)) Ngrid_box_larger = j
     if( abs(rae(j) -biggerboxradius) <  abs(rae(Ngrid_biggerbox) -biggerboxradius)) Ngrid_biggerbox = j
  end do

  Ngrid_box_larger=Ngrid_biggerbox
  

  !!  if(abs(rae(Ngrid_box) -boxradius)>1.0e-8) STOP "the grid from pseudo should pass by rcov but it does not"
  
  iovrsmpl=1 + (Ngrid-1)/ngrid_fit
  !! iovrsmpl=1 + (Ngrid-1)/Ngrid_biggerbox 
  
  
  Ngrid =  (iovrsmpl*(ngrid_fit-1)+1)
  !! Ngrid =  (iovrsmpl*(Ngrid_biggerbox-1)+1)
  
  allocate( rgrid   (Ngrid ))
  allocate( rw  (Ngrid ))
  allocate( rd   (Ngrid ))
  allocate( nonloc   (Ngrid ))
  allocate( rgrid_ab(Ngrid ))
  allocate( aepot(Ngrid ))
  allocate( aepot_cent(Ngrid ))
  
  allocate( staepot(Ngrid ))


  !! IN ATOM :    r(i) = a*(exp(b*(i-1))-1)
  b=log( (rae(201)-rae(101))/(rae(101)-rae(1)))/100.0_8
  a= rae(201)/( exp(b*(201-1))-1)
  
  !! IN PAWPATCH : more points
  b=b/iovrsmpl            
  

  fourpi=16.d0*atan(1.d0)
  do igrid=1,Ngrid
     !! rgrid( igrid ) =  a*(exp(b*(igrid-1))-1)
     rgrid( igrid ) =  a*(exp(b*(igrid-1))    )          !!  uniform logarithmic grid, tending to a*(exp(b*(igrid-1))-1)
                                                       !! otherwise the logarithmic step
                                                       !! is not constant and gets too big
                                                       !! going to r=0
     !! rgrid_ab(igrid) = (rgrid(igrid)+a) * b     NOT used in schro
     rw(igrid)=b*rgrid(igrid)
     rd(igrid)=1.d0/rw(igrid)
     rw(igrid)=rw(igrid)*fourpi*rgrid(igrid)**2
  end do
  
  Ngrid_box= 1+iovrsmpl*( Ngrid_box-1)
  Ngrid_biggerbox=1+iovrsmpl*( Ngrid_biggerbox-1)
  Ngrid_box_larger = Ngrid_biggerbox

  boxradius = rgrid(Ngrid_box)
  

  !!  if(abs(rgrid(Ngrid_box) -boxradius)>1.0e-8) STOP "the finer grid should still pass by rcov but it does not"

    
  allocate(yp (ngrid_fit))
  allocate(ypp(ngrid_fit))
  allocate(w(3*ngrid_fit))

 
  allocate( aepot_p (Ngrid ))
  allocate( aepot_pp(Ngrid ))
    
  a1=0.0D0
  b1=0.0D0
  aN=0.0D0
  bN=0.0D0
  
  isx=0
  call splift(rae,   atom_potential_fit,yp, ypp, ngrid_fit,w, ierr,isx, a1,b1,aN,bN )
 

  write(plotfile, '(a,i0,a)') 'ae.pot.spline.',nconfpaw ,'.plt'
  open(unit=37,file=trim(plotfile),status='unknown')
  do j=1,ngrid_fit
     write(37,*)  rae(j),atom_potential_fit (j), yp(j),ypp (j)
  end do
  close(37)  

  call splint(rae,atom_potential_fit, ypp, ngrid_fit, &
       rgrid, aepot, aepot_p, aepot_pp, Ngrid, ierr)

  if ( trim(pawstatom)/='') then
     isx=0
     call splift(rae,statom_potential,yp, ypp, ngrid_fit,w, ierr,isx, a1,b1,aN,bN )
     call splint(rae,statom_potential, ypp, ngrid_fit, &
          rgrid, staepot, aepot_p, aepot_pp, Ngrid, ierr)
     open(unit=37,file="stae.pot.interp",status='unknown')
     do j=1,Ngrid
        write(37,*)  rgrid(j), staepot(j)
     end do
     close(37)  
  endif
  
  deallocate(aepot_p  )
  deallocate(aepot_pp )
  deallocate(yp)
  deallocate(ypp)
  deallocate(w)

  write(plotfile, '(a,i0,a)') 'ae.pot.interp.',nconfpaw ,'.plt'
  open(unit=37,file=trim(plotfile),status='unknown')
  do j=1,Ngrid
     write(37,*)  rgrid(j), aepot(j)
  end do
  close(37)  
  

  EMAX =  (Nsol*2*3.1415/boxradius)**2 *1.5 !! to have a margin

  aepot = aepot +15.0D0
  if ( trim(pawstatom)/='') staepot = staepot +15.0D0

  if(ispp=='r') then
     do iout=1,2
        uout=outunits(iout)
        write(uout,*)  "WARNING : relativistic calculation but pawpatch will use non relativistic solution " 
        write(uout,*)  "WARNING : A future correction of this problem would be using the Koelling-Hammon solution " 
        write(uout,*)  "WARNING : to the scalar  dirac equation J. Phys. C: Solid State Phys., Vol. 10. 1977 " 
     end do
  endif

  allocate( psigrid(Ngrid  , Nsol ))
  allocate( Egrid  (Nsol         ))

  allocate( psigrid_pseudo(Ngrid  , Nsol ))
  allocate( psipsigrid_pseudo(Ngrid  , Nsol ))
  allocate( Egrid_pseudo  (Nsol         ))

  allocate( psigrid_bigger(Ngrid  , Nsol ))

  allocate(psi_initial_copy( Ngrid )) 
  allocate(psi_initial( Ngrid )) 
  allocate(dumpsi_p(Ngrid))
  allocate(expo(ng))
  allocate( PAWpatch_matrix(Npaw,Npaw ))

  nonloc=0.0_8
  if ( trim(pawstatom)/='') then

     aepot_cent= 0.5_8*staepot + (pawstL*(pawstL+1) )/rgrid/rgrid/2   !! energies were in Rydeberg
     
     call schro( dum_energy ,rgrid , &
          aepot_cent  ,nonloc, psi_initial_copy(1)   , Ngrid ,&
          pawstN  , pawstL  ,   znuc )
     
     ! dum_energy=-2000.0D0 
     ! call difnrl(aepot_cent ,psi_initial_copy(1) , dumpsi_p  ,&
     !      Ngrid, a,b, rgrid,rgrid_ab, pawstN+pawstL , pawstL , znuc,dum_energy, EMAX )
     
     
     write(38,*)  "initial wave L=",pawstL ," N= " , pawstN,   " E= ", dum_energy
     write(6,* )    , "initial wave L=",pawstL ," N= " , pawstN,   " E= ", dum_energy
     
     write(filename,'(a,I0)')'psi_initial_L_',pawstL
     open(unit=22,file=trim(filename))
     do igrid=1, Ngrid
        write(22,'(2(f20.10,1x))') rgrid(igrid),psi_initial_copy(igrid )  
     enddo
     close(unit=22)
     
     do igrid=1, Ngrid
        psi_initial_copy(igrid)=psi_initial_copy(igrid)*rgrid(igrid)**pawstP 
     enddo
  endif
  

  write(6,*)  "CALCULATING PAWpatch corrections for l from 0 to " , npawl-1
  write(6,*)  "------------------------------------------------------ " 
  write(38,*) "CALCULATING PAWpatch corrections for l from 0 to " , npawl-1
  write(38,*) "------------------------------------------------------ "

  do LPaw=0, npawl-1     
     psi_initial=psi_initial_copy

     write(6,*)   "==============================================================="
     write(6,*)   "========== now CALCULATING PAWpatch correction  for l  =  " ,  LPaw
     write(6,*)   "==============================================================="
     write(38,*)   "==============================================================="
     write(38,*)   "========= now CALCULATING PAWpatch correction  for l  =  " ,  LPaw
     write(38,*)   "==============================================================="

     aepot_cent=aepot*0.5_8 + (Lpaw*(Lpaw+1) )/rgrid/rgrid/2 

     psigrid=0.0D0
 
     write(6,*)     "now calculating " , 10, " AE function for the eigenvalues CHECK with  Ngrid_box_larger  for LPaw=", LPaw 
     write(38,*)  "now calculating " , 10, " AE function for the eigenvalues CHECK with  Ngrid_box_larger  for LPaw=", LPaw 

     do n=1, 10
        !! call difnrl(aepot_cent ,psigrid(1, n ) , dumpsi_p ,Ngrid_box, a,b, rgrid,rgrid_ab,n+Lpaw ,Lpaw , znuc,Egrid(n) )
        call schro( Egrid(n),rgrid , &
             aepot_cent  ,nonloc,  psigrid(1, n )  , Ngrid_box_larger ,&
             n+LPaw , Lpaw  ,   znuc )                     
        write(6 ,*)  " n = ", n, " Egrid(n)  "  , Egrid(n),  "  lpaw " , Lpaw
        write(38,*)  "schro ae  n = ", n, " Egrid(n)  for CHECK with   " , Egrid(n),  "  lpaw " , Lpaw
     enddo
    
     write(6,*)     "now calculating " , NSol, " function of the AE basis for LPaw=", LPaw 
     write(38,*)  "now calculating " , NSol, " function of the AE basis for LPaw=", LPaw 

     do n=1, Nsol
        !! call difnrl(aepot_cent ,psigrid(1, n ) , dumpsi_p ,Ngrid_box, a,b, rgrid,rgrid_ab,n+Lpaw ,Lpaw , znuc,Egrid(n) )
        call schro( Egrid(n),rgrid , &
             aepot_cent  ,nonloc,  psigrid(1, n )  , Ngrid_box ,&
             n+LPaw , Lpaw  ,   znuc )                     
        write(38,*)  "schro ae  n = ", n, " Egrid(n) " , Egrid(n),  "  lpaw " , Lpaw
        if( psigrid(Ngrid_box-1, n)<0) then
           do igrid=1, ngrid
              psigrid(igrid, n)=-psigrid(igrid,n)
           end do
        endif
     enddo
 
     if( dump_functions) then
        write(plotfile, '(a,i0,a)') 'ae.wfs.L=',LPaw,'.plt'
        open(unit=22,file= trim(plotfile) )
        do igrid=1, Ngrid
           write(22,'(200(f20.10,1x))') rgrid(igrid), (psigrid(igrid,j ), j=1,Nsol)   
        enddo
        close(unit=22)
     endif

     !! rimetto energia in Hartree
     !! Egrid=Egrid/2


     Egrid_pseudo(:)= Egrid(:)  !! to fit these energies and find the dual
     !! Egrid different from zero, with psp_modifier=0,  activates the fit of psigrid
     !! energy by energy and the calculation of paw stuff

     psigrid_pseudo=psigrid
 

     call paw_generator(znuc,zion,lmx,lpmx,lmax,hsep, gpot, &
          rloc, r_l, &                         !! rloc=alpz=alpl    r_l=alps
          ng-1 ,noccmax ,noccmx,   expo,  psi,aeval, occup ,  &
          Nsol, Lpaw , Ngrid, Ngrid_box,Egrid_pseudo,  rgrid , rw,rd, psigrid_pseudo ,&
          Npaw, PAWpatch_matrix,  psipsigrid_pseudo, rcov, rprb, rcore,zcore, Ngrid_box_larger)
     
     if( dump_functions) then
        write(plotfile, '(a,i0,a)') 'ptildes.L=',LPaw,'.plt'
        open(unit=22,file=trim(plotfile))
        do igrid=1, Ngrid
           write(22,'(200(f20.10,1x))') rgrid(igrid), (psigrid_pseudo(igrid,j ), j=1,Npaw)   
        enddo
        close(unit=22)
        
        write(plotfile, '(a,i0,a)') 'psitildes.L=',LPaw,'.plt'
        open(unit=22,file=trim(plotfile))
        do igrid=1, Ngrid
           write(22,'(200(f20.10,1x))') rgrid(igrid), (psipsigrid_pseudo(igrid,j ), j=1,Npaw)   
        enddo
        close(unit=22)
     endif
     
     
     real_start=-1
     !! write(38,*)  ">>  Comparaison Egrid  Egrid_pseudo "
     !! write(6 ,*)  ">>  Comparaison betwenn first 5 Egrid and   Egrid_pseudo "
     do n=1, Nsol
        if((Egrid(n)+0.1).ge.Egrid_pseudo(1)) then
           if(real_start==-1) then
               real_start = n
            end if
            !! write(38,*) Egrid(n), Egrid_pseudo(n -real_start+1)
            !! if (n.le.5) write(6 ,*) Egrid(n), Egrid_pseudo(n -real_start+1)
        else
           !! write(38,*) Egrid(n)
           !! if (n.le.5) write(6 ,*) Egrid(n)
        endif
     enddo


     !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     !!  In the following   AE are calculated on a bigger grid to 
     !! avoid compression of the  missing bound states
     !!
     
     psigrid_bigger=0.0D0
     do n=1, Nsol
        
        call schro( Egrid(n),rgrid , &
             aepot_cent  ,nonloc,  psigrid_bigger(1, n )  , Ngrid_biggerbox ,&
             n +LPaw , Lpaw  ,   znuc )

        ! call difnrl(aepot_cent ,psigrid_bigger(1, n ) , dumpsi_p ,Ngrid_biggerbox, a,b, &
        !      rgrid,rgrid_ab,n+Lpaw ,Lpaw , znuc,), EMAX )
     enddo
     
     
     if(iproc.eq.0 .and. dump_functions) then
        write(plotfile, '(a,i0,a)') 'wf_bigger.L=',LPaw,'.plt'
        open(unit=22,file=trim(plotfile))
        do igrid=1, Ngrid
           write(22,'(200(f20.10,1x))') rgrid(igrid), (psigrid_bigger(igrid,j ), j=1,Nsol)   
        enddo
        close(unit=22)
     endif
     

     if ( trim(pawstatom)/='') then
        write(6,*) "routine pawpatch  , PROJECT  initial wf*r**pawstP on pseudos "
        write(38,*)  "routine pawpatch  , PROJECT  initial wf*r**pawstP on pseudos "
        call find_pfproj_4tail( Nsol,Npaw,Ngrid,  Ngrid_box,Ngrid_biggerbox, rgrid, psi_initial, psigrid, real_start, &
             psigrid_pseudo, psipsigrid_pseudo,  &
             psigrid_bigger,dump_functions) 
        

        !! if(iproc.eq.0 .and. dump_functions.eq.1) then 
        write(plotfile, '(a,i0,a)') 'projres.L=',LPaw,'.plt'
        open(unit=22,file=trim(plotfile))
        do igrid=1, Ngrid_biggerbox
           !! la funzione proiettata e in ultima colonna
           write(22,'(200(f20.10,1x))') rgrid(igrid),  psi_initial(igrid), psigrid(igrid,1), psigrid(igrid,2)
        enddo
        close(unit=22)
     endif


     do iout=1,2
        uout=outunits(iout)
        write(uout,*) "Initial Projected wf at igrid=100,500,2000"
        write(uout,*) psigrid(100,2)
        write(uout,*) psigrid(500,2)
        write(uout,*) psigrid(2000,2)
     end do
        
     write(plotfile, '(a,i0,a)') 'pawdata.L=',LPaw,'.dat'
     open(unit=22,file=trim(plotfile))
          
     write(6,*)  "now writing the dual functions and PAWpatch for ", Npaw," duals "
     write(22,'(A)') "now writing the dual functions and PAWpatch  "
     write(22,'(I4,1x,I6,1x,I4)') Npaw, Ngrid_box, Lpaw
     do igrid=1, Ngrid_box
        write(22,'(E20.14)') rgrid(igrid)
     enddo
     
     write(22,'(A)') " "
     do n=1, Npaw
        do igrid=1, Ngrid_box
           write(22,'(E20.14)') psigrid_pseudo(igrid, n)
        enddo
        write(22,'(A)') " "
     enddo
     

     write(6,*) "PawPatch Matrix"
     do n=1, Npaw
           write(6,*) PawPatch_matrix(n,:)
     enddo
     
     do n=1, Npaw
        do j=1, Npaw
           write(22,'(E20.14)') PawPatch_matrix(j,n)
        enddo
     enddo
     
     
     if(pawstatom/=" " ) then
        write(6,*)  "now writing the dual functions and PAWpatch for ", 1," duals for initial wave "
        write(22,'(A)') "now writing the dual functions and PAWpatch  "
        write(22,'(I4,1x,I6,1x,I4)') 1, Ngrid_biggerbox, -Lpaw-1
        do igrid=1, Ngrid_biggerbox
           write(22,'(E20.14)') rgrid(igrid)
        enddo
        
        write(22,'(A)') " "
        do n=1, 1
           do igrid=1, Ngrid_biggerbox
              write(22,'(E20.14)') psigrid(igrid,2)
           enddo
           write(22,'(A)') " "
        enddo
        
        do n=1, 1
           do j=1, 1
              write(22,'(E20.14)') 0.0_8
           enddo
        enddo
     endif

     close(unit=22)
     
  enddo
  deallocate(Egrid)
  deallocate(psigrid)
  deallocate(Egrid_pseudo)
  deallocate(psigrid_pseudo)
  deallocate(psi_initial)
  deallocate(dumpsi_p)
  deallocate(atom_potential_fit)
  deallocate(statom_potential)
  deallocate(staepot)
  deallocate(expo)
  deallocate(PAWpatch_matrix)
  
  stop
END subroutine pawpatch



subroutine find_pfproj_4tail( Nsol,Npaw, Ngrid,Ngrid_box,Ngrid_biggerbox,&
     rgrid, psi1s, psigrid, real_start,&
     ptilde, psitilde, &
     psigrid_bigger, dump_functions )
  implicit none
  integer, parameter :: gp=kind(1.0d0) 
  !Arguments
  integer, intent(in) ::  Nsol,Npaw,Ngrid,Ngrid_box,Ngrid_biggerbox,real_start
  real(gp), intent(inout) :: psi1s(Ngrid), rgrid(Ngrid)
  real(gp), intent(inout) :: psigrid(Ngrid,Nsol),psigrid_bigger(Ngrid,Nsol)
  real(gp), intent(inout) :: psitilde(Ngrid,Nsol), ptilde(Ngrid,Nsol)
  !! real(gp) , intent(out) :: coeffs_out(Npaw)
  logical :: dump_functions
  !Local variables
  real(gp)   :: coeffs_out(Npaw)
  real(gp) :: dumgrid(Ngrid),  dumgrid2(Ngrid)
  integer :: i,k
  real(gp)  :: coeffs(Nsol), ratio, x

  !! check
  do i=1, Nsol-real_start+1

     dumgrid =psigrid_bigger(:,i)*psigrid_bigger(:,i)
     call integrate(dumgrid, dumgrid2, rgrid, Ngrid_biggerbox)
     if( abs(dumgrid2(Ngrid_biggerbox)-1.0_gp).gt.1.0D-5) Then
        write(6,*)  "  norm(psigrid_bigger) != 1 in find_pfproj_4tail"
        STOP" program stopped, problem in find_pfproj_4tail" 
     endif


     dumgrid =psigrid(:,i)*psigrid(:,i)
     call integrate(dumgrid, dumgrid2, rgrid, Ngrid_box)
     if( abs(dumgrid2(Ngrid_box)-1.0_gp).gt.1.0D-5) Then
        write(6,*)  "  norm(psigrid) != 1 in find_pfproj_4tail", dumgrid2(Ngrid_box)-1.0_gp
        STOP" program stopped, problem in find_pfproj_4tail" 
     endif

     dumgrid =psigrid_bigger(:,1)*psigrid(:,i)
     call integrate(dumgrid, dumgrid2, rgrid, Ngrid_box)
     !! print * , " >>>>>>>> " , i, "  " , dumgrid2(Ngrid_box)

  end do

  do i=1, real_start-1

     do k=1, Ngrid
        dumgrid(k)=psigrid_bigger(k,i)*psi1s(k)
     enddo

     call integrate(dumgrid, dumgrid2, rgrid, Ngrid_biggerbox)


     coeffs(i)=dumgrid2(Ngrid_biggerbox)
     write(6,*)  "From initial wave  subtract psigrid_bigger(",i,") with coeff ", coeffs(i)  
     write(38,*) "From initial wave  subtract AE psigrid_bigger(",i, ") with coeff ", coeffs(i)  
     do k=1, Ngrid
        psi1s(k) = psi1s(k) -coeffs(i)*psigrid_bigger(k,i)
     enddo
  end do

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
  end do
 

  coeffs_out(:)=coeffs(real_start:real_start+Npaw-1)


  do i=1, Nsol-real_start+1
     ratio = psigrid( Ngrid_box-10,i+real_start-1)/psitilde( Ngrid_box-10,i)
     !! print *, "psigrid ", i , "would require q correction factor ", ratio 
     !!$ psitilde(:,i)=ratio *psitilde(:,i)
  enddo

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


  psigrid(:,1)=dumgrid
  psigrid(:,2)=dumgrid2


  return
END SUBROUTINE find_pfproj_4tail
