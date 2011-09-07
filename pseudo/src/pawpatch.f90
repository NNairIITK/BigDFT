subroutine pawpatch(energ,verbose,maxdim,pp,penal,&
     noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,pol,nsmx,&
     no,lo,so,ev,crcov,dcrcov,ddcrcov,norb,&
     occup,aeval,chrg,dhrg,ehrg,res,wght,&
     wfnode,psir0,wghtp0,&
     rcov,rprb,rcore,zcore,znuc,zion,rloc,gpot,r_l,hsep,&
     vh,xp,rmt,rmtg,ud,nint,ng_fit,ngmx,psi,&
     avgl1,avgl2,avgl3,ortprj,litprj,igrad,rr_fit,rw,rd, &
     iproc,nproc,wghtconf,wghtexci,wghtsoft,wghtrad,wghthij,&
     nhgrid,hgridmin,hgridmax, nhpow,ampl,crmult,frmult,&
     excitAE,ntime,iter,itertot,penref,time,ngrid_fit,&
     nconfpaw, npawl, nchannelspaw, ispp, pawstatom , pawstN, pawstL, pawstP )
  
  implicit none
  !! implicit real*8 (a-h,o-z)
  logical avgl1,avgl2,avgl3,ortprj,litprj,igrad
  real(8) pp(maxdim),so(norb),ev(norb),crcov(norb),&
       dcrcov(norb),ddcrcov(norb),occup(noccmx,lmx,nsmx),aeval(noccmx,lmx,nsmx),&
       chrg(noccmx,lmx,nsmx),dhrg(noccmx,lmx,nsmx),&
       ehrg(noccmx,lmx,nsmx),res(noccmx,lmx,nsmx),&
       wght(noccmx,lmx,nsmx,8),&
       wfnode(noccmx,lmx,nsmx,3),&
       gpot(*),r_l(*),hsep(6,lpmx,nsmx),&
       vh(*),xp(*),rmt(*),rmtg(*),ud(*),psi(*),&
       rr_fit(*),rw(*),rd(*),pen_cont(7),&
       exverbose(4*nproc),time(3), penal, psir0,wghtp0,rcov,&
       rprb,rcore,zcore,znuc,zion,rloc,&
       wghtexci,wghtsoft,wghtrad,wghthij,&
       hgridmin,hgridmax, ampl,crmult,frmult,&
       excitAE,iter,itertot,penref, wghtconf
  
  integer no(norb),lo(norb),nconfpaw, npawl, nchannelspaw , maxdim,&
       noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx,nint,ng_fit,ngmx,iproc,&
       nproc,nhgrid,nhpow,ntime, norb, ngrid_fit, j,   pawstN, pawstL, pawstP
  
  character(len=1) ispp
  logical:: energ, verbose, wpen, pol
  character(len=30) :: plotfile
  real(8), pointer :: atom_potential_fit(:)
  real(8), pointer :: statom_potential(:)
  real(8) rdum
  
  integer Nsolm, Npaw, ng
  integer Ngrid, Ngrid_box, Ngrid_biggerbox, iovrsmpl
  real(8) boxradius, biggerboxradius, a,b
  real(8), pointer :: rgrid(:), yp(:), ypp(:), w(:), aepot(:), aepot_p(:), aepot_pp(:), &
       rgrid_ab(:), aepot_cent(:)
  real(8) a1,b1,an,bn
  integer ierr, isx
  integer LPaw, n, Nsol
  integer igrid
  real(8), pointer :: psi_initial(:), dumpsi_p(:)
  real(8) dum_energy
  character(1000) filename
  character(len=125) :: pawstatom
  real(gp) , pointer ::  psigrid(:,:), Egrid(:)
  real(gp) , pointer ::  psigrid_pseudo(:,:), Egrid_pseudo(:)
  real(gp) , pointer ::  expo(:)
  real(gp), pointer::PAWpatch(:,:)
  logical dumpfunctions
  
  dumpfunctions= .true.



  include 'mpif.h'
  

  if (nspin/=1) then
     stop " pawpatch can be built only in the nspin=1 case   " 
  endif

  if(iproc/=0) return
  print *, "IN PAWPATCH   ngrid_fit=", ngrid_fit 
  allocate(atom_potential_fit(ngrid_fit))
  write(plotfile, '(a,i0,a)') 'ae.pot.conf.',nconfpaw ,'.plt'
  open(unit=37,file=trim(plotfile),status='unknown')
  do j=1,ngrid_fit
     read(37, *)   rdum, atom_potential_fit(j)
     if( abs(rdum-rr_fit(j))>1.0e-6) then
        STOP "rgrid not corresponding "
     end if
  end do
  close(37)
  
  allocate(statom_potential(ngrid_fit))
  if ( trim(pawstatom)/='') then
     open(unit=37,file=trim(pawstatom),status='unknown')
     do j=1,ngrid_fit
        read(37, *)   rdum, statom_potential(j)
        if( abs(rdum-rr_fit(j))>1.0e-6) then
           STOP "rgrid in statom  not corresponding  with the one in ae.pot.conf "
        end if
     end do
     close(37)
  endif

  ng  = 30
  !! noccmax = 5 
  print *, "  NOCCMAX as set by pseudo is  ", noccmax
  !! lmax=3
  print *,"LMAX as set by pseudo is ", lmax 
  
  Nsol=200
  
  Npaw= nchannelspaw
  
  Ngrid=20000                      !! roughly the number of point of the 
  !! oversampled grid
  
  boxradius=rcov                  !! this should be found in the original grid
  biggerboxradius = 1.5_8 * rcov   !! this is an approximative goal
  
  Ngrid_box=1
  Ngrid_biggerbox=1
  
  do j=1,ngrid_fit
     if( abs(rr_fit(j) -boxradius) <  abs(rr_fit(Ngrid_box) -boxradius)) Ngrid_box = j
     if( abs(rr_fit(j) -biggerboxradius) <  abs(rr_fit(Ngrid_biggerbox) -biggerboxradius)) Ngrid_biggerbox = j
  end do
  if(abs(rr_fit(Ngrid_box) -boxradius)>1.0e-8) STOP "the grid from pseudo should pass by rcov but it does not"
  
  !! iovrsmpl=1 + (Ngrid-1)/ngrid_fit
  iovrsmpl=1 + (Ngrid-1)/Ngrid_biggerbox 
  
  
  !! Ngrid =  (iovrsmpl*(ngrid_fit-1)+1)
  Ngrid =  (iovrsmpl*(Ngrid_biggerbox-1)+1)
  
  allocate( rgrid   (Ngrid ))
  allocate( rgrid_ab(Ngrid ))
  allocate( aepot(Ngrid ))
  allocate( aepot_cent(Ngrid ))
  
  allocate( staepot(Ngrid ))


  !! IN ATOM :    r(i) = a*(exp(b*(i-1))-1)
  b=log( (rr_fit(201)-rr_fit(101))/(rr_fit(101)-rr_fit(1)))/100.0_8
  a= rr_fit(201)/( exp(b*(201-1))-1)
  
  !! IN PAWPATCH : more points
  b=b/iovrsmpl            
  
  do igrid=1,Ngrid
     rgrid( igrid ) =  a*(exp(b*(igrid-1))-1)
     rgrid_ab(igrid) = (rgrid(igrid)+a) * b
  end do
  
  Ngrid_box= 1+iovrsmpl*( Ngrid_box-1)
  Ngrid_biggerbox=1+iovrsmpl*( Ngrid_biggerbox-1)

  if(abs(rgrid(Ngrid_box) -boxradius)>1.0e-8) STOP "the finer grid should still pass by rcov but it does not"

    
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
  call splift(rr_fit,   atom_potential_fit,yp, ypp, ngrid_fit,w, ierr,isx, a1,b1,aN,bN )
 

  write(plotfile, '(a,i0,a)') 'ae.pot.spline.',nconfpaw ,'.plt'
  open(unit=37,file=trim(plotfile),status='unknown')
  do j=1,ngrid_fit
     write(37,*)  rr_fit(j),atom_potential_fit (j), yp(j),ypp (j)
  end do
  close(37)  

  call splint(rr_fit,atom_potential_fit, ypp, ngrid_fit, &
       rgrid, aepot, aepot_p, aepot_pp, Ngrid, ierr)

  if ( trim(pawstatom)/='') then
     isx=0
     call splift(rr_fit,statom_potential,yp, ypp, ngrid_fit,w, ierr,isx, a1,b1,aN,bN )
     call splint(rr_fit,statom_potential, ypp, ngrid_fit, &
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
  
  aepot = aepot +15.0D0-1000.0D0
  if ( trim(pawstatom)/='') staepot = staepot +15.0D0-1000.0D0
  
  if(ispp=='r') then
     print *, " WARNING : relativistic calculation but pawpatch will use non relativistic solution " 
     print *, " WARNING : A future correction of this problem would be using the Koelling-Hammon solution " 
     print *, " WARNING : to the scalar  dirac equation J. Phys. C: Solid State Phys., Vol. 10. 1977 " 
  endif

  allocate( psigrid(Ngrid  , Nsol ))
  allocate( Egrid  (Nsol         ))

  allocate( psigrid_pseudo(Ngrid  , Nsol ))
  allocate( Egrid_pseudo  (Nsol         ))

  allocate(psi_initial( Ngrid )) 
  allocate(dumpsi_p(Ngrid))
  deallocate(expo(ng))
  allocate( PAWpatch(Npaw,Npaw ))

  do LPaw=0, npawl-1     
     if ( trim(pawstatom)/='') then
        dum_energy=-2000.0D0
        aepot_cent=staepot + (Lpaw*(Lpaw+1) )/rgrid/rgrid   !! energies are in Rydeberg
        call difnrl(aepot_cent ,psi_initial(1) , dumpsi_p ,&
             Ngrid, a,b, rgrid,rgrid_ab,n+Lpaw ,Lpaw , znuc,dum_energy )
        print *, "L=", LPaw,  " E= ", dum_energy
        write(filename,'(a,I0)')'psi_initial_L_',LPaw
        open(unit=22,file=trim(filename))
        do igrid=1, Ngrid
           write(22,'(4(f20.10,1x))') rgrid(igrid), (psi_initial(igrid,j ), j=1,3)   
        enddo
        close(unit=22)
     endif

     aepot_cent=aepot + (Lpaw*(Lpaw+1) )/rgrid/rgrid 

     psigrid=0.0D0
     do n=1, Nsol
        call difnrl(aepot_cent ,psigrid(1, n ) , dumpsi_p ,Ngrid_box, a,b, rgrid,rgrid_ab,n+Lpaw ,Lpaw , znuc,Egrid(n) )
     enddo
     if( dump_functions) then
        write(plotfile, '(a,i0,a)') 'ae.wfs.L=',LPaw,'.plt'
        open(unit=22,file= trim(plotfile) )
        do igrid=1, Ngrid
           write(22,'(200(f20.10,1x))') rgrid(igrid), (psigrid(igrid,j ), j=1,Nsol)   
        enddo
        close(unit=22)
     endif

     
     Egrid_pseudo(:)= Egrid(:)  !! to fit these energies and find the dual
     !! Egrid different from zero, with psp_modifier=0,  activates the fit of psigrid
     !! energy by energy and the calculation of paw stuff
     print *, "copy psigrid_pseudo "
     psigrid_pseudo=psigrid
!!$  print *, "  chiamo   abs_generator_modified con psipsigrid_pseudo, Ngrid_biggerbox  ", Ngrid_biggerbox
!!$ psipsigrid_pseudo(1,1)=Ngrid_biggerbox


     call paw_generator( znuc , zion   , lpmx,    hsep, gpot, &
          rloc, r_l, &                         !! rloc=alpz=alpl    r_l=alps
          ng-1 ,noccmax ,  expo,  psi,aeval, occup ,  &
          Nsol, abs_final_L , Ngrid, Ngrid_box,Egrid_pseudo,  rgrid , psigrid_pseudo ,&
          Npaw, PAWpatch,  psipsigrid_pseudo)
     

     if( dump_functions.eq.1) then
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
     print *, " Comparaison Egrid  Egrid_pseudo "
     do iocc=1, Nsol
        if((Egrid(iocc)+0.1).ge.Egrid_pseudo(1)) then
           if(real_start==-1) then
              real_start = iocc
              print *, Egrid(iocc)
           else
              print *, Egrid(iocc), Egrid_pseudo(iocc -real_start+1)
           endif
        endif
     enddo

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
  deallocate(PAWpatch)

  stop
END subroutine pawpatch

