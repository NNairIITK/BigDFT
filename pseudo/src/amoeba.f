

c     minimizes the penalty function (subroutine) using
c     a simplex downhill method (amoeba)

c     cwh: if amoeba doesn't improve within ntrymax iterations
c     amoeba also returns

      SUBROUTINE AMOEBA(P,Y,ndim,FTOL,ITER,itmax,namoeb,
     :     noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,pol,nsmx,
     :     no,lo,so,ev,crcov,dcrcov,ddcrcov,norb,
     :     occup,aeval,chrg,dhrg,ehrg,res,wght,
     :     wfnode,psir0,wghtp0,
     :     rcov,rprb,rcore,zcore,znuc,zion,rloc,gpot,r_l,hsep,
     :     vh,xp,rmt,rmtg,
     :     ud,nint,ng,ngmx,psi,avgl1,avgl2,avgl3,
     :     ortprj,litprj,igrad,rr,rw,rd,
c          the following line differs from pseudo2.2
     :     iproc,nproc,wghtconf,wghtexci,wghtsoft,wghtrad,wghthij,
     :     nhgrid, hgridmin,hgridmax, nhpow,ampl,crmult,frmult,
     :     ntrymax,excitAE,ntime,itertot,energ,verbose,time)



      implicit real*8 (a-h,o-z)
      logical  avgl1,avgl2,avgl3,ortprj,litprj,igrad,
     :         lexit,lnext,verbose,pol
      PARAMETER (NMAX=50,ALPHA=1.0d0,BETA=0.5d0,GAMMA=2.0d0)
c      parameter ( trymax = 200 ) is now a dummy argument read from FITPAR

      DIMENSION P(ndim,ndim+1),Y(ndim+1),PR(NMAX),PRR(NMAX),PBAR(NMAX)
      dimension no(*),lo(*),so(*),ev(*),crcov(*),dcrcov(*),ddcrcov(*),
     :     occup(*),aeval(*),chrg(*),dhrg(*),ehrg(*),res(*),wght(*),
     :     gpot(*),r_l(*),hsep(*),vh(*),xp(*),rmt(*),rmtg(*),
     :     ud(*),psi(*),wfnode(*),rr(*),rw(*),rd(*),time(3)
      include 'mpif.h'

c      print*,'entered amoeba with nfit=',ndim,itmax,FTOL,ntrymax
      if (ndim.gt.nmax-1) stop 'nmax'
      if (ndim.lt.1) then
         write(6,*)'entered amoeba with nfit=',ndim
         write(6,*)'no fit!'
         return
      endif
       
      MPTS=NDIM+1
      ITER=0
      ntrycount = 0
      iloold = 1

c     ================== HERE begins the loop for the simplex ======================

c     first, get the currently highest, second highest and lowest vertices of the simplex

 1    ILO=1
      IF(Y(1).GT.Y(2))THEN
         IHI=1
         INHI=2
      ELSE
         IHI=2
         INHI=1
      ENDIF
      DO 11 I=1,MPTS
         IF(Y(I).LT.Y(ILO)) ILO=I
         IF(Y(I).GT.Y(IHI))THEN
            INHI=IHI
            IHI=I
         ELSE IF(Y(I).GT.Y(INHI))THEN
            IF(I.NE.IHI) INHI=I
         ENDIF


 11   CONTINUE

c     Then do various checks for exit conditions

c     cwh
      if (ilo .ne. iloold) ntrycount = 0
c     write(6,*)'debug: ilo, ntrycount,ntrymax',ilo, ntrycount,ntrymax
      iloold = ilo

c     RTOL=min(Y(IHI),(y(ihi)-y(ilo))/y(ilo)**2)
      RTOL=min(Y(IHI),(y(ihi)-y(ilo))/y(ilo))
      IF(RTOL.LT.FTOL) then
c     Check
         write(6,*) 'Converged:',rtol,'<',ftol
c        write(6,*) 'values at the edges of the simplex:'
c        write(6,'(40e15.8)') y
         do i=1,ndim
            if (y(i).lt.y(ilo)) write(6,*) 'WARNING ilo not lowest'
         enddo
      endif
c      if (mod(iter,100).eq.0) write(6,*) 'iter=',iter
c
      IF(ITER.gt.ITMAX) then
c        write(6,*) 'WARNING: NO CONVERGENCE IN AMOEBA'
         write(6,*) '--- Simplex done, max no of iterations reached.'
c         ftol=10.d0*ftol
c         write(6,*) 'FTOL SET TO ',FTOL
c         write(6,*) 'values at the edges of the simplex:'
c         write(6,'(40e15.8)') y
c         write(6,*)'ilo:',ilo
         do i=1,ndim
            if (y(i).lt.y(ilo)) write(6,*) 'WARNING ilo not lowest'
         enddo
      end if

      IF(ntrycount.ge.ntrymax) then
         write(6,*)'no improvement during',ntrycount,
     :             'iterations, tolerance is exceeded.'
c        write(6,*) 'WARNING: NO IMPROVEMENT IN AMOEBA FOR THE LAST',
c    :        NTRYCOUNT,'STEPS'
      ENDIF

c     every 10th step, call a barrier and check for exit/next files
      if (mod(iter,10).eq.0) then
c        a barrier for synchronisation cannot be wrong here
         if (nproc > 1) call MPI_BARRIER(MPI_COMM_WORLD,ierr)  
         INQUIRE ( FILE = 'NEXT', EXIST = lnext )
         INQUIRE ( FILE = 'EXIT', EXIST = lexit )
         if (lnext) then 
            if (nproc > 1) call MPI_BARRIER(MPI_COMM_WORLD,ierr)  
            write(6,*) 'The file NEXT is present, aborting this fit' 
            if(iproc.eq.0)then
               open(99,file='NEXT')
               close(99, status='DELETE')
            end if
         endif
         if (lexit) then
c           set namoeb to something smaller than iiter
c           so that the main programs cycle loop exits 
            namoeb=-namoeb
            if (nproc > 1) call MPI_BARRIER(MPI_COMM_WORLD,ierr)  
            write(6,*) 'The file EXIT is present, aborting fits' 
            if(iproc.eq.0)then
               open(99,file='EXIT')
               close(99, status='DELETE')
            end if
         endif
      endif


c     Several exit conditions to leave amoeba, but always the same actions:
c             pack the lowest vertex in the current psppar
c             call penalty to make sure the latest call to gatom was indeed
c             for the lowest vertex and therefore the packed psppar
c             call penalty rather than gatom to give some more information
c             about the excitation energies and the softness when leaving
c             the simplex cycle

      IF( (ITER.gt.ITMAX) .or. lexit .or. lnext .or.
     :    (ntrycount.ge.ntrymax) .or. (RTOL.LT.FTOL) )then



c          another call to penalty with the verbose flag yields further information

           write(6,*)'______________'
           write(6,*)'leaving amoeba'
           write(6,*)
           write(6,*)
           verbose=.true.
c          write(6,*)'verbose=',verbose,' inquire penalty information'
           

           call penalty(energ,verbose,ndim,p(1,ilo),y(ilo),
     :        noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,pol,nsmx,
     :        no,lo,so,ev,crcov,dcrcov,ddcrcov,norb,
     :        occup,aeval,chrg,dhrg,ehrg,res,wght,
     :        wfnode,psir0,wghtp0,
     :        rcov,rprb,rcore,zcore,znuc,zion,rloc,gpot,r_l,hsep,
     :        vh,xp,rmt,rmtg,ud,nint,ng,ngmx,psi,
     :        avgl1,avgl2,avgl3,ortprj,litprj,igrad,rr,rw,rd,
c             the following line differs from pseudo2.2
     :        iproc,nproc,wghtconf,wghtexci,wghtsoft,wghtrad,wghthij,
     :        nhgrid, hgridmin,hgridmax, nhpow,ampl,crmult,frmult,
     :        excitAE,ntime,iter,itertot,Y(ilo),time)


c             unpack variables from p(ilo) into psppar variables
              verbose=.false.
              call  ppack (verbose,rloc,gpot,hsep,r_l,p(1,ilo),
     :              lpx,lpmx,nspin,pol,nsmx,ndim,ndim,'unpack',avgl1,
     :              avgl2,avgl3,ortprj,litprj,
     :              rcore,zcore,znuc,zion)
c             verbose=.false.
         RETURN
      endif


C     CWH
      ntrycount = ntrycount +1
      if(mod(ntrycount,max(1,ntrymax/5))==0)
     :   write(6,*)'no improvement during',ntrycount,
     :             'iterations, tolerance is',ntrymax

c     The simplex is monitored with the penalty contributions in penalty.f
c        write(6,*) 'iter (amoeba,gatom):',iter,ntime,
c    :        ' y(ilo):',y(ilo)
c???         CALL FLUSH(6)




c     ================== HERE are the actual simplex downhill moves ======================

      ITER=ITER+1
c     get the simplex centre
      DO 12 J=1,NDIM
         PBAR(J)=0.d0
 12   CONTINUE
      DO 14 I=1,MPTS
         IF(I.NE.IHI)THEN
            DO 13 J=1,NDIM
               PBAR(J)=PBAR(J)+P(j,i)
 13         CONTINUE
         ENDIF
 14   CONTINUE

c     reflect the highest vertex at the centre -> PR
      DO 15 J=1,NDIM
         PBAR(J)=PBAR(J)/NDIM
         PR(J)=(1.d0+ALPHA)*PBAR(J)-ALPHA*P(j,IHI)
 15   CONTINUE
      call penalty(energ,verbose,ndim,pr,ypr,
     :     noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,pol,nsmx,
     :     no,lo,so,ev,crcov,dcrcov,ddcrcov,norb,
     :     occup,aeval,chrg,dhrg,ehrg,res,wght,
     :     wfnode,psir0,wghtp0,
     :     rcov,rprb,rcore,zcore,znuc,zion,rloc,gpot,r_l,hsep,
     :     vh,xp,rmt,rmtg,
     :     ud,nint,ng,ngmx,psi,avgl1,avgl2,avgl3,
     :     ortprj,litprj,igrad,rr,rw,rd,
c          the following line differs from pseudo2.2
     :     iproc,nproc,wghtconf,wghtexci,wghtsoft,wghtrad,wghthij,
     :     nhgrid, hgridmin,hgridmax, nhpow,ampl,crmult,frmult,
     :     excitAE,ntime,iter,itertot,Y(ilo),time)
      IF(YPR.LE.Y(ILO))THEN
c           new lowest, can we go further in that direction?

c       test a bigger step for the reflecting 
        DO 16 J=1,NDIM
          PRR(J)=GAMMA*PR(J)+(1.d0-GAMMA)*PBAR(J)
16      CONTINUE
        call penalty(energ,verbose,ndim,prr,yprr,
     :       noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,pol,nsmx,
     :       no,lo,so,ev,crcov,dcrcov,ddcrcov,norb,
     :       occup,aeval,chrg,dhrg,ehrg,res,wght,
     :       wfnode,psir0,wghtp0,
     :       rcov,rprb,rcore,zcore,znuc,zion,rloc,gpot,r_l,hsep,
     :       vh,xp,rmt,rmtg,
     :       ud,nint,ng,ngmx,psi,avgl1,avgl2,avgl3,
     :       ortprj,litprj,igrad,rr,rw,rd,
c            the following line differs from pseudo2.2
     :       iproc,nproc,wghtconf,wghtexci,wghtsoft,wghtrad,wghthij,
     :       nhgrid, hgridmin,hgridmax, nhpow,ampl,crmult,frmult,
     :       excitAE,ntime,iter,itertot,Y(ilo),time)
        IF(YPRR.LT.Y(ILO))THEN
c         this got even better, so accept PRR, discarding the highest vertex 
          DO 17 J=1,NDIM
            P(j,IHI)=PRR(J)
17        CONTINUE
c          if (mod(iter,10).eq.0) write(6,'(a,i5,a,e15.7,a,i2,e15.7)')
c     :         'iter',iter,' found',YPRR,' reject',ihi,Y(IHI)
          Y(IHI)=YPRR
        ELSE
c         PR is the best we have, so accept it, discarding the highest vertex 
          DO 18 J=1,NDIM
            P(j,IHI)=PR(J)
18        CONTINUE
c          if (mod(iter,10).eq.0) write(6,'(a,i5,a,e15.7,a,i2,e15.7)')
c     :         'iter',iter,' found',YPR,' reject',ihi,Y(IHI)
          Y(IHI)=YPR
        ENDIF
      ELSE IF(YPR.GE.Y(INHI))THEN
c       the reflected vertex is not lower than the second highest vertex
        IF(YPR.LT.Y(IHI))THEN
c         if reflecting improved the highest vertex, update it
          DO 19 J=1,NDIM
            P(j,IHI)=PR(J)
19        CONTINUE
c          if (mod(iter,10).eq.0)
c     :         write(6,'(a,i5,a,e15.7,2(a,i2,e15.7))')
c     :         'iter',iter,' found',YPR,' reject',ihi,Y(IHI),
c     :         ' best:',ilo,Y(Ilo)
          Y(IHI)=YPR
        ENDIF
c       try to contract the highest vertex to the centre
        DO 21 J=1,NDIM
          PRR(J)=BETA*P(j,IHI)+(1.d0-BETA)*PBAR(J)
21      CONTINUE
      call penalty(energ,verbose,ndim,prr,yprr,
     :     noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,pol,nsmx,
     :     no,lo,so,ev,crcov,dcrcov,ddcrcov,norb,
     :     occup,aeval,chrg,dhrg,ehrg,res,wght,
     :     wfnode,psir0,wghtp0,
     :     rcov,rprb,rcore,zcore,znuc,zion,rloc,gpot,r_l,hsep,
     :     vh,xp,rmt,rmtg,
     :     ud,nint,ng,ngmx,psi,avgl1,avgl2,avgl3,
     :     ortprj,litprj,igrad,rr,rw,rd,
c          the following line differs from pseudo2.2
     :     iproc,nproc,wghtconf,wghtexci,wghtsoft,wghtrad,wghthij,
     :     nhgrid, hgridmin,hgridmax, nhpow,ampl,crmult,frmult,
     :     excitAE,ntime,iter,itertot,Y(ilo),time)
        IF(YPRR.LT.Y(IHI))THEN
c         if contracting improved the highest vertex, save it
          DO 22 J=1,NDIM
            P(j,IHI)=PRR(J)
22        CONTINUE
c          if (mod(iter,10).eq.0)
c     :         write(6,'(a,i5,a,e15.7,2(a,i2,e15.7))')
c     :         'iter',iter,' found',YPRR,' reject',ihi,Y(IHI),
c     :         ' best:',ilo,Y(Ilo)
          Y(IHI)=YPRR
        ELSE
c       all moves of the highest vertex failed to improof the penalty,
c       so we shrink the entire simplex towards the lowest vertex
          DO 24 I=1,MPTS
            IF(I.NE.ILO)THEN
              DO 23 J=1,NDIM
                PR(J)=0.5d0*(P(j,I)+P(j,ILO))
                P(j,I)=PR(J)
23            CONTINUE
      call penalty(energ,verbose,ndim,pr,ypr,
     :     noccmax,noccmx,lmax,lmx,lpx,lpmx,lcx,nspin,pol,nsmx,
     :     no,lo,so,ev,crcov,dcrcov,ddcrcov,norb,
     :     occup,aeval,chrg,dhrg,ehrg,res,wght,
     :     wfnode,psir0,wghtp0,
     :     rcov,rprb,rcore,zcore,znuc,zion,rloc,gpot,r_l,hsep,
     :     vh,xp,rmt,rmtg,
     :     ud,nint,ng,ngmx,psi,avgl1,avgl2,avgl3,
     :     ortprj,litprj,igrad,rr,rw,rd,
c          the following line differs from pseudo2.2
     :     iproc,nproc,wghtconf,wghtexci,wghtsoft,wghtrad,wghthij,
     :     nhgrid, hgridmin,hgridmax, nhpow,ampl,crmult,frmult,
     :     excitAE,ntime,iter,itertot,Y(ilo),time)
            ENDIF
24        CONTINUE
        ENDIF
      ELSE
        DO 25 J=1,NDIM
          P(j,IHI)=PR(J)
25      CONTINUE
c          if (mod(iter,10).eq.0)
c     :         write(6,'(a,i5,a,e15.7,2(a,i2,e15.7))')
c     :         'iter',iter,' found',YPR,' reject',ihi,Y(IHI),
c     :         ' best:',ilo,Y(Ilo)
        Y(IHI)=YPR
      ENDIF
      GO TO 1
      END



