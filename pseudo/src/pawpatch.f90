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
      real(8), pointer :: rgrid

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
      iovrsmpl=1 + (Ngrid-1)/ngrid_fit
      Ngrid =  (iovrsmpl*(ngrid_fit-1)+1)
      
      allocate( rgrid(Ngrid ))

      !! IN ATOM :    r(i) = a*(exp(b*(i-1))-1)
      b=log( (rr_fit(201)-rr_fit(101))/(rr_fit(101)-rr_fit(1)))/100.0_8
      a= rr_fit(201)/( exp(b*(201-1))-1)

      !! IN PAWPATCH : more points
      b=b/iovrsmpl            
      Ngrid_box= 1+iovrsmpl*( Ngrid_box-1)
      Ngrid_biggerbox=1+iovrsmpl*( Ngrid_biggerbox-1)
      do i=1,Ngrid
            rgrid( i ) =  a*(exp(b*(i-1))-1)
            rgrid_ab(i)= (rgrid(i)+a)*b
      end do
      if(abs(rgrid(Ngrid_box) -boxradius)>1.0e-8) STOP "the finer grid should still pass by rcov but it does not"
      
      call splift(rr_fit,  atom_potential_fit





      deallocate(atom_potential_fit)
      stop
    END subroutine pawpatch
    








! c     set res() =-1 for orbitals with zero charge and
! c      wght(nocc,l+1,ispin,5) set to zero
! c     this avoids unneccessay computations of res() in the
! c     routine resid()

! c     when we print out detailed info with verbose,
! c     we actually do not want to miss any residues

!       do iorb=1,norb
!          nocc=no(iorb)
!          l=lo(iorb)
!          ispin=1
!          if (so(iorb).lt.0) ispin=2
!          if ( (wght(nocc,l+1,ispin,5).eq.0.0d0) .and.
!      :        (occup(nocc,l+1,ispin).lt.1.0d-8) .and.
!      :        ( .not. verbose) )  then
!             res(nocc,l+1,ispin) = -1.d0
!          else
!             res(nocc,l+1,ispin) = 0.d0
!          endif
!       enddo
! c  unpack variables
! c      print*,'penalty: maxdim=',maxdim
! c     convention for spin treatment with libXC using the pol variable
!             nspol=1
!             if(pol)nspol=2
!       call  ppack (verbose,rloc,gpot,hsep,r_l,pp(1),
!      :     lpx,lpmx,nspin,pol,nsmx,maxdim,maxdim,'unpack',
!      :     avgl1,avgl2,avgl3,ortprj,litprj,
!      :     rcore,zcore,znuc,zion)

!       call cpu_time(t)
!       time(1)=time(1)-t
!       call gatom(nspol,energ,verbose,
!      :     noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx,
!      :     occup,aeval,chrg,dhrg,ehrg,res,wght,wfnode,psir0,
!      :     rcov,rprb,rcore,zcore,znuc,zion,rloc,gpot,r_l,hsep,
!      :     vh,xp,rmt,rmtg,ud,nint,ng,ngmx,psi,igrad,
!      :     rr,rw,rd,ntime,itertot,etotal)
!       call cpu_time(t)
!       time(1)=time(1)+t
!       penal=0d0
! c     diff for dcharg and echarge is calc. in (%)
!       do iorb=1,norb
!          nocc=no(iorb)
!          l=lo(iorb)
!          ispin=1
!          if (so(iorb).lt.0) ispin=2
!          penal=penal + 
!      :        ((aeval(nocc,l+1,ispin)-ev(iorb))
!      :        *wght(nocc,l+1,ispin,1))**2 +
!      :        ((chrg(nocc,l+1,ispin)-crcov(iorb))
!      :        *wght(nocc,l+1,ispin,2))**2 +
!      :        (100.d0*(1.d0-dhrg(nocc,l+1,ispin)/dcrcov(iorb))
!      :        *wght(nocc,l+1,ispin,3))**2 +
!      :        (100.d0*(1.d0-ehrg(nocc,l+1,ispin)/ddcrcov(iorb))
!      :        *wght(nocc,l+1,ispin,4))**2  

!               aa=wfnode(nocc,l+1,ispin,1)
!               if(.not. (aa>0d0.or.aa<=0d0))then
!                      write(6,*)'WARNING, wfnode=NaN for iorb',iorb
!                      cycle
!               end if
!          penal=penal + 
!      :        (res(nocc,l+1,ispin)*wght(nocc,l+1,ispin,5))**2 +
!      :        (wfnode(nocc,l+1,ispin,1)*wght(nocc,l+1,ispin,6))**2 +
!      :        (wfnode(nocc,l+1,ispin,2)*wght(nocc,l+1,ispin,7))**2 +
!      :        (wfnode(nocc,l+1,ispin,3)*wght(nocc,l+1,ispin,8))**2  


! c             write(66,*)'DEBUG: iorb, nocc,l,ispin',iorb, nocc,l,ispin
! c             write(66,*)'evals',aeval(nocc,l+1,ispin),ev(iorb)
! c             write(66,*)'DEBUG: node', wfnode(nocc,l+1,ispin,1)
! c             write(66,*)'DEBUG, all wghts',wght(nocc,l+1,ispin,1:8)
! c             write(66,*)'DEBUG, penalty sum: iorb, penal',iorb, penal
! c             write(66,*)
!       enddo
! c     for now, keep this contribution here, no special output
!       penal=penal + (psir0*wghtp0)**2

! c     add weight for narrow radii
! c     use: error in norm conversation of projectors
! c          is approximately a power law of the radius
! c
! c         1-|projector|  ~ 1d-6 (hgrid/r)**12

! c     write(16,*)'DEBUG: sqpenal without rad',penal
!       pen_r=0d0
!       do l=1,lpx+1
!           pen_r=pen_r+(hgridmax/r_l(l))**24
!       end do
! c     just use the same convention for rloc
!       pen_r=pen_r+(hgridmax/rloc)**24
! c     prefactors
!       pen_r=pen_r*1d-12*wghtrad**2
! c     write(16,*)'DEBUG: sqpenal with rad',penal

! c     add an empirical term to penalize all hsep larger than ~10
! c     also calculate kij and penalize values larger than ~ 2,
! c     otherwise we might favor big splittings.
! c     Use a power law such that h~10 and k~2 contribute each about
! c     1% of wghthij to the penalty, but larger values are strictly
! c     penalized.
!       pen_h=0d0
!       pen_k=0d0
!       do l=1,lpx+1
!         do i=1,6
!            if(hsep(i,l,1)==0)cycle
!            pen_h= pen_h+(hsep(i,l,1)*0.1d0)**12
!            if(hsep(i,l,2)==0)cycle
!            pen_h= pen_h+(hsep(i,l,2)*0.1d0)**12
!            tk=2d0*(hsep(i,l,1)-hsep(i,l,2))/(2*l-1)
!            if(tk==0)cycle
! c          we may want to try different weights on kij later
!            pen_k= pen_k+(tk*0.5d0)**12
!         end do
!       end do
!       if(verbose)then
!         write(6,*)
!         write(6,*)'Empirical penalty terms from psppar'
!         write(6,*)'___________________________________'
!         write(6,*)
!         write(6,*)'radii   ',sqrt(pen_r)
!         write(6,*)'hij     ',sqrt(pen_h)*wghthij*1d-2
!         write(6,*)'kij     ',sqrt(pen_k)*wghthij*1d-2
!         write(6,*)
     
!       end if
           
!       pen_cont(7)=pen_r+(pen_h+pen_k)*1d-4*wghthij**2



! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! c---------------- parallel computation of excitation energies and softness ----------------
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! c     first, make sure that excitation energies of value zero are ignored
!       if(excitAE==0d0)wghtexci=0d0
      
! c     if the weight for softness is nonzero,
! c     calculate the change in kinetic energy when
! c     transforming PSI into a debauchie wavelet basis

!       if(wghtsoft > 0d0 .and. nhgrid>0 )then ! .and. mod(iter,nskip)==0)
!          call ekin_wvlt(verbose,iproc,nproc,ng,ngmx,
!      :   noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,nsmx,
!      :   nhgrid, hgridmin,hgridmax, nhpow,ampl,crmult,
!      :        frmult, xp,psi,occup,ekin_pen,time)
! c        important: Here we scale the radii for the coarse
! c                   and fine grid with the parameter rcov.
!       else
!          ekin_pen=0d0
!       end if

!       if(nproc==1)then
! c        Serial case: One configuration, no excitation energies
!          pen_cont(4)=wghtconf**2*penal
!          pen_cont(6)=wghtsoft**2*ekin_pen
!          penal=sqrt(penal+sum(pen_cont(6:7)))
! c        write(16,*)'DEBUG:ekin_pen,penal',ekin_pen,penal,wghtsoft
!       else 
! c        parallel case
! c        Set up penalty contributions to be summed over all processes
!          call cpu_time(t)
!          time(2)=time(2)-t

! c        compute excitation energy for this configuration. 
! c        for this, get total energy of configuration 1.
!          if(iproc==0) eref=etotal
!          call MPI_BCAST(eref,1,MPI_DOUBLE_PRECISION,0,
!      :                         MPI_COMM_WORLD,ierr)
!          if(ierr.ne.0)write(*,*)'MPI_BCAST ierr',ierr
!          excit=etotal-eref


! c        Sum up penalty contributions from this MPI process
!          pen_cont(1)=wghtconf**2 *penal
!          pen_cont(2)=wghtexci**2 *(excit-excitAE)**2
!          pen_cont(3)=wghtsoft**2 *ekin_pen

! c        write(6,*)'DEBUG lines for penalty over processes'
! c        write(6,*)'wghtconf,penal',wghtconf ,penal
! c        write(6,*)'wghtexci,(excit-excitAE)**2',
! c    :       wghtexci,(excit-excitAE)**2
! c        write(6,*)'wghtsoft,ekin_pen',wghtsoft,ekin_pen
! c        write(6,*)pen_cont(1:3),'products'

! c        sum up over all processes, that is over all configurations
! c        and over the orbitals and hgrid samples of the wavelet transform
!          call MPI_ALLREDUCE(pen_cont(1:3),pen_cont(4:6),3,
!      :         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!          if(ierr.ne.0)write(*,*)'MPI_ALLREDUCE ierr',ierr
!          call cpu_time(t)
!          time(2)=time(2)+t
!          penal=sqrt(sum(pen_cont(4:7)))

!       end if


! c     If the penalty is below the reference - usually the currently
! c     lowest vertex of the simplex -  write out the major components:


!       if(penal<penref)then
!           if(nproc==1)then
!              write(6,'(2i8,4(1x,e20.10))')
!      :           iter, ntime, penal, sqrt(pen_cont(6)),
!      :           sqrt(pen_cont(7)),sqrt(pen_cont(4))
!           else
!              pen_cont(4:7)=sqrt(pen_cont(4:7))
!              write(6, '(2i8,7(1x,e20.10))') 
!      :            iter, ntime, penal,pen_cont(6),pen_cont(7),
!      :            pen_cont(4:5),sqrt(sum(pen_cont(1:2)))
!           end if
! c         write the vertex to the dumpfile without transforming  psppar
! c         either only the lowest or its history 
! c         backspace(99)
!           if(iproc==0)then
!             write(99,'(6e11.3,t65,a)')
!      :         penal,rloc,gpot(1),gpot(2),gpot(3),gpot(4),
!      :         'penalty, rloc, gpot'
!                do l=0,lpx
!                   write(99,'(f7.3,t8,6e11.3,t76,a)') r_l(l+1),
!      :                 (hsep(i,l+1,1),i=1,6),'r_l(),hsep(up)'
!                   if (l.gt.1-nspol .and. nspin.eq.2)
!      :                 write(99,'(t8,6e11.3,t76,a)')
!      :                 (hsep(i,l+1,2),i=1,6),'      hsep(dn)'
!                enddo
!             write(99,*)  
!           end if
!       end if





! c
!       if(nproc>1.and.verbose)then
! c         write out the excitation energies
! c         get excitation energies from all processes  
! c         using mpiallreduce may be a clumsy way of doing this
!           exverbose=0d0
!           exverbose(2*iproc+1)=excit
!           exverbose(2*iproc+2)=sqrt(pen_cont(2))
! c         note: Write out energies in the convention they have been read
!           call MPI_ALLREDUCE(exverbose(1),
!      :            exverbose(2*nproc+1),2*nproc,
!      :            MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!           write(6,*)
!           write(6,*)'Excitation energies'
!           write(6,*)'___________________'
!           write(6,*)
!           write(6,*)'Configuration,    dE=Etot-Etot1,'//
!      :              '    (dE-dE_AE)*weight'
!           do i=1,nproc
!              write(6,'(10x,i4,3e20.12)')
!      :                 i-1,exverbose(2*nproc+2*i-1),
!      :                 exverbose(2*nproc+2*i)
!           end do
!           write(6,*)
!       end if
!       ierr=0
!       if(nproc>1) call MPI_BARRIER(MPI_COMM_WORLD,ierr)  
!       if(ierr.ne.0)write(*,*)'MPI_BARRIER ierr',ierr

!       return
!       end

