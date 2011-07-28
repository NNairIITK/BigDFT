           subroutine pawpatch(energ,verbose,maxdim,pp,penal,&
           noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,pol,nsmx,&
           no,lo,so,ev,crcov,dcrcov,ddcrcov,norb,&
           occup,aeval,chrg,dhrg,ehrg,res,wght,&
           wfnode,psir0,wghtp0,&
           rcov,rprb,rcore,zcore,znuc,zion,rloc,gpot,r_l,hsep,&
           vh,xp,rmt,rmtg,ud,nint,ng,ngmx,psi,&
           avgl1,avgl2,avgl3,ortprj,litprj,igrad,rr_fit,rw,rd, &
           iproc,nproc,wghtconf,wghtexci,wghtsoft,wghtrad,wghthij,&
           nhgrid,hgridmin,hgridmax, nhpow,ampl,crmult,frmult,&
           excitAE,ntime,iter,itertot,penref,time,ngrid_fit,&
           nconfpaw, npawl, nchannelspaw )






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
           noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx,nint,ng,ngmx,iproc,&
           nproc,nhgrid,nhpow,ntime, norb, ngrid_fit, j

      

      logical:: energ, verbose, wpen, pol
      character(len=30) :: plotfile
      real(8), pointer :: atom_potential_fit(:)
      real(8) rdum

      integer Nsolm Npaw, ng, noccmax, lmax
      integer Ngrid, Ngrid_box, Ngrid_biggerbox, iovrsmpl
      real(8) boxradius, biggerboxradius, a,b
      real(8), pointer :: rgrid(:), yp(:), ypp(:), w(:), aepot(:), aepot_p(:), aepot_pp(:)
      real(8) a1,b1,an,bn
      integer ierr, isx
      integer LPaw, n
      real(8), pointer :: psi_initial(:,:), dumpsi_p(:)
      real(8) dum_energy
      character(1000) filename
      include 'mpif.h'
      
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

      ng  = 30
      noccmax = 5 
      lmax=3
      
      Nsol=200
      
      Npaw= nchannelspaw

      Ngrid=20000                      !! roughly the number of point of the 
                                       !! oversampled grid
  
      boxradius=rcov                   !! this should be found in the original grid
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
      
      allocate( rgrid(Ngrid ))
      allocate( aepot(Ngrid ))

      !! IN ATOM :    r(i) = a*(exp(b*(i-1))-1)
      b=log( (rr_fit(201)-rr_fit(101))/(rr_fit(101)-rr_fit(1)))/100.0_8
      a= rr_fit(201)/( exp(b*(201-1))-1)

      !! IN PAWPATCH : more points
      b=b/iovrsmpl            

      do i=1,Ngrid
            rgrid( i ) =  a*(exp(b*(i-1))-1)
            rgrid_ab(i)= (rgrid(i)+a)*b
      end do
     
    
      Ngrid_box= 1+iovrsmpl*( Ngrid_box-1)
      Ngrid_biggerbox=1+iovrsmpl*( Ngrid_biggerbox-1)
      if(abs(rgrid(Ngrid_box) -boxradius)>1.0e-8) STOP "the finer grid should still pass by rcov but it does not"
 




      allocate(yp (ngrid_fit))
      allocate(ypp(ngrid_fit))
      allocate(w(3*ngrid_fit))

      a1=0.0D0
      b1=0.0D0
      aN=0.0D0
      bN=0.0D0
      
      call splift(rr_fit,   atom_potential_fit,yp, ypp, ngrid_fit,w, ierr,isx, a1,b1,aN,bN )

      allocate( aepot_p (Ngrid ))
      allocate( aepot_pp(Ngrid ))

      call splint(rr_fit,atom_potential_fit, ypp, ngrid_fit, &
           rgrid, aepot, aepot_p, aepot_pp, Ngrid, ierr)

      deallocate(aepot_p  )
      deallocate(aepot_pp )
      deallocate(yp)
      deallocate(ypp)
      deallocate(w)


      
      aepot = aepot +15.0D0-1000.0D0
  
      do LPaw=0, npawl-1     

         allocate(psi_initial( Ngrid, 3))  !! up to thre different n : 1s  2s  3s.. 2p  3p  4p... ( if it makes sense )
         allocate(dumpsi_p(Ngrid))
         do n=1,3
            call difnrl(aepot ,psi_initial(1,n) , dumpsi_p ,&
                 Ngrid, a,b, rgrid,rgrid_ab,n+Lpaw ,Lpaw , znuc,dum_energy )
            print *, "L=", LPaw, " n=", n, " E= ", dum_energy
         end do
         deallocate(dumpsi_p)
         
         write(filename,'(a,a0)')'numerov_initial_L_',LPaw
         open(unit=22,file=trim(filename))
         do igrid=1, Ngrid
            write(22,'(4(f20.10,1x))') rgrid(igrid), (psi_initial(igrid,j ), j=1,3)   
         enddo
         close(unit=22)
      enddo



      deallocate(atom_potential_fit)
      stop
    END subroutine pawpatch
    
