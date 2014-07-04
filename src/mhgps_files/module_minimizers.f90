module module_minimizers

contains

!> @file
!!  Routines for Stefan's new minimization method
!! @author
!!    Copyright (C) 2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
!subroutine geopt(nat,wpos,etot,fout,fnrmtol,count,count_sd,displr)
subroutine minimizer_sbfgs(runObj_,outsIO,verbosity,ncount_bigdft,fail)
!call_bigdft has to be run once on runObj_ and outs !before calling this routine
!sbfgs will return to caller the energies and coordinates used/obtained from the last accepted iteration step
   use module_base
   use module_types
   use module_interfaces
   use yaml_output
   use module_global_variables
   use module_energyandforces
   implicit none
   !parameter
   integer, intent(in)                    :: verbosity
   type(run_objects), intent(inout)       :: runObj_
   type(DFT_global_output), intent(inout) :: outsIO
   integer, intent(inout)                 :: ncount_bigdft
   logical, intent(out)                   :: fail
   !local variables
   character(len=*), parameter :: subname='sbfgs'
   integer :: info !< variables containing state codes
   integer :: nhistx !< maximum history length
   integer :: nhist  !< actual history length
   integer :: ndim   !< dimension of significant subspace
   integer :: nit    !< maximum number of iterations
   integer :: nat    !< number of atoms
   integer :: istat,iall
   integer :: it,i,iat,l,j,idim,jdim,ihist,icheck !<counter variables
   integer :: itswitch
!   type(DFT_global_output) :: outs
   logical :: debug !< set .true. for debug output to fort.100
   logical :: steep !< steepest descent flag
   real(gp) :: displr !< (non-physical) integrated path length,
                      !< includes displacements from rejctions
                      !< (denoted as dsplr in geopt.mon)
   real(gp) :: displp !< (physical) integrated path length,
                      !< includes NO displacements from rejections
                      !< (denoted as dsplp in geopt.mon)
   real(gp) :: etot 
   real(gp) :: fnrm
   real(gp) :: fmax
   real(gp) :: fluct
   real(gp) :: betax !< initial step size (gets not changed)
   real(gp) :: beta  !< current step size
   real(gp) :: cosangle
   real(gp) :: ts
   real(gp) :: tt
   real(gp) :: maxd !< maximum displacement of single atom
   real(gp) :: scl
   real(gp) :: etotold
   real(gp) :: detot
   real(gp) :: etotp
   real(gp) :: dt
   real(gp) :: cutoffRatio !< if fraction of eigenvalues of overlapmatrix of
                           !< displacements is smaller that cutoffRatio,
                           !< then those displacements are regarded
                           !< as linear dependent and are not taken into account
                           !< for building the significant subspace
   real(gp) :: maxrise !< energy ist not allowed to rise more than maxrise in single step
   real(gp) :: steepthresh !< if fnrm is larger that steepthresh, steepest descent is used
   real(gp) :: trustr !< a single atoms is not allowed to be dsiplaced more than by trustr
   real(gp), allocatable, dimension(:,:,:) :: rxyz
   real(gp), allocatable, dimension(:,:,:) :: fxyz
   real(gp), allocatable, dimension(:,:,:) :: ff
   real(gp), allocatable, dimension(:,:,:) :: rr
   real(gp), allocatable, dimension(:,:,:) :: rrr
   real(gp), allocatable, dimension(:,:,:) :: fff
   real(gp), allocatable, dimension(:,:)   :: aa
   real(gp), allocatable, dimension(:,:)   :: dd
   real(gp), allocatable, dimension(:)     :: eval
   real(gp), allocatable, dimension(:)     :: res
   real(gp), allocatable, dimension(:)     :: scpr
   real(gp), allocatable, dimension(:)     :: rnorm
   character(len=4)                        :: fn4
   character(len=40)                       :: comment
   character(len=9)                        :: cdmy9_1
   character(len=9)                        :: cdmy9_2
   character(len=9)                        :: cdmy9_3
   character(len=8)                        :: cdmy8


   !set parameters
   nit=runObj%inputs%ncount_cluster_x
   nat=runObj%atoms%astruct%nat
   betax=runObj%inputs%betax
   nhistx=runObj%inputs%nhistx
   maxrise=runObj%inputs%maxrise
   cutoffRatio=runObj%inputs%cutoffratio
   steepthresh=runObj%inputs%steepthresh
   trustr=runObj%inputs%trustr
   if (iproc==0.and.verbosity > 0) then
      call yaml_open_map('Geometry parameters')
         call yaml_map('Geometry Method','GEOPT_SBFGS')
         call yaml_map('nhistx',nhistx)
         call yaml_map('betax', betax,fmt='(1pe21.14)')
         call yaml_map('maxrise', maxrise,fmt='(1pe21.14)')
         call yaml_map('cutoffRatio', cutoffRatio,fmt='(1pe21.14)')
         call yaml_map('steepthresh', steepthresh,fmt='(1pe21.14)')
         call yaml_map('trustr', trustr,fmt='(1pe21.14)')
      call yaml_close_map()
   end if

   !init varaibles
   debug=.false.
   fail=.true.
   displr=0.0_gp
   displp=0.0_gp
   fluct=0.0_gp
   icheck=0
   detot=0.0_gp
   itswitch=0
   ndim=0
   nhist=0
   beta=betax
   maxd=0.0_gp

   ! allocate arrays
   rxyz = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='rxyz')
   fxyz = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='fxyz')
   aa = f_malloc((/ nhistx, nhistx /),id='aa')
   eval = f_malloc(nhistx,id='eval')
   res = f_malloc(nhistx,id='res')
   rnorm = f_malloc(nhistx,id='rnorm')
   ff = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='ff')
   rr = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='rr')
   dd = f_malloc((/ 3, nat /),id='dd')
   fff = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='fff')
   rrr = f_malloc((/ 1.to.3, 1.to.nat, 0.to.nhistx /),id='rrr')
   scpr = f_malloc(nhistx,id='scpr')
!   call init_global_output(outs, runObj%atoms%astruct%nat)

   !copy outs_datatype
!   call copy_global_output(outsIO,outs)



!!!!!!   call energyandforces(nat,rxyz(1,1,0),fxyz(1,1,0),etot)
!!  not necessary, call_bigdft allready called outside
!   call call_bigdft(runObj,outs,nproc,iproc,infocode)
!   ncount_bigdft=ncount_bigdft+1

!! copy to internal variables
   call vcopy(3*runObj%atoms%astruct%nat, runObj%atoms%astruct%rxyz(1,1), 1,rxyz(1,1,0), 1)
   call vcopy(3*outs%fdim, outs%fxyz(1,1), 1, fxyz(1,1,0), 1)
   etot=outs%energy

   call fnrmandforcemax(fxyz(1,1,0),fnrm,fmax,nat)
   fnrm=sqrt(fnrm)
   if (fmax < 3.d-1) call updatefluctsum(outs%fnoise,fluct)

   etotold=etot
   etotp=etot

   if (iproc==0.and.verbosity > 0) then
       !avoid space for leading sign (numbers are positive, anyway)
       write(cdmy8,'(es8.1)')abs(maxd)
       write(cdmy9_1,'(es9.2)')abs(displr)
       write(cdmy9_2,'(es9.2)')abs(displp)
       write(cdmy9_3,'(es9.2)')abs(beta)

       write(16,'(i5,1x,i5,2x,a10,2x,1es21.14,2x,es9.2,es11.3,3es10.2,2x,a6,a8,xa4,i3.3,xa5,a7,2(xa6,a8))') &
       ncount_bigdft,0,'GEOPT_SBFGS',etotp,detot,fmax,fnrm,fluct*runObj%inputs%frac_fluct,fluct, &
       'beta=',trim(adjustl(cdmy9_3)),'dim=',ndim,'maxd=',trim(adjustl(cdmy8)),'dsplr=',trim(adjustl(cdmy9_1)),'dsplp=',trim(adjustl(cdmy9_2))
   endif

   do it=1,nit!start main loop
!  do it=1,nit-1!start main loop (nit-1 if first bigdft call is NOT done outside, but inside this subroutine)
      if (debug.and.iproc==0) write(100,*) 'it:',it,etot,fnrm,itswitch
      nhist=nhist+1

      if (fnrm.gt.steepthresh .or. it.le.itswitch ) then
         ndim=0
         steep=.true.
         if (it.gt.itswitch) itswitch=it+nhistx
         if (debug.and.iproc==0) write(100,*) "STEEP"
      else
         steep=.false.
      endif

      ! make space in the history list
      if (nhist.gt.nhistx) then
         nhist=nhistx
         do ihist=0,nhist-1
            do iat=1,nat
               do l=1,3
                  rxyz(l,iat,ihist)=rxyz(l,iat,ihist+1)
                  fxyz(l,iat,ihist)=fxyz(l,iat,ihist+1)
               enddo
            enddo
         enddo
      endif
   
      ! decompose gradient
500 continue
      do iat=1,nat
         do l=1,3
            dd(l,iat)=-fxyz(l,iat,nhist-1)
         enddo
      enddo
      do i=1,ndim
         scpr(i)=0.0_gp
         do iat=1,nat
            do l=1,3
               scpr(i)=scpr(i)-fxyz(l,iat,nhist-1)*rrr(l,iat,i)
            enddo
         enddo
         do iat=1,nat
            do l=1,3
               dd(l,iat)=dd(l,iat)-scpr(i)*rrr(l,iat,i)
            enddo
         enddo
      enddo
   
      ts=0.0_gp
      do iat=1,nat
         do l=1,3
            ts=ts+dd(l,iat)**2
         enddo
      enddo
         
      if (debug.and.iproc==0) write(100,*) 'beta=',beta
      do iat=1,nat
         do l=1,3
            dd(l,iat)=dd(l,iat)*beta
         enddo
      enddo
   
      do i=1,ndim
         !eval(i) is corrected by possible error (Weinstein criterion)
         tt=scpr(i)/sqrt(eval(i)**2+res(i)**2)
         if (debug.and.iproc==0) write(100,'(a,i3,3(1x,e10.3))') 'i,tt,eval,res '&
                    ,i,1.0_gp/sqrt(eval(i)**2+res(i)**2),eval(i),res(i)
            do iat=1,nat
               do l=1,3
                  dd(l,iat)=dd(l,iat)+tt*rrr(l,iat,i)
               enddo
            enddo
      enddo
   
      tt=0.0_gp
      dt=0.0_gp
      maxd=-huge(1.0_gp)
      do iat=1,nat
         dt=dd(1,iat)**2+dd(2,iat)**2+dd(3,iat)**2
         tt=tt+dt
         maxd=max(maxd,dt)
      enddo
      tt=sqrt(tt)
      maxd=sqrt(maxd)
   
      !trust radius approach: avoids too large steps due to large forces
      !only used when in steepest decent mode
      if(maxd>trustr .and. steep)then
         if(debug.and.iproc==0)write(100,'(a,xes24.17,xi0)')'step too large',maxd,it
         if(iproc==0)write(16,'(a,2(xes9.2))')'WARNING GEOPT_SBFGS: step too large: maxd, trustradius ',maxd,trustr
         scl=0.50_gp*trustr/maxd
         dd=dd*scl
         tt=tt*scl
         maxd=maxd*scl
      endif
      displr=displr+tt
   
      !update positions
      do iat=1,nat
         rxyz(1,iat,nhist)=rxyz(1,iat,nhist-1)-dd(1,iat)
         rxyz(2,iat,nhist)=rxyz(2,iat,nhist-1)-dd(2,iat)
         rxyz(3,iat,nhist)=rxyz(3,iat,nhist-1)-dd(3,iat)
      enddo
   
      runObj%inputs%inputPsiId=1
      call energyandforces(nat,runObj%atoms%astruct%cell_dim,rxyz(1,1,nhist),fxyz(1,1,nhist),etotp)
!      call energyandforces(nat,rxyz(1,1,nhist),fxyz(1,1,nhist),etotp)
!      call vcopy(3 * runObj%atoms%astruct%nat, rxyz(1,1,nhist), 1,runObj%atoms%astruct%rxyz(1,1), 1)
!      runObj%inputs%inputPsiId=1
!      call call_bigdft(runObj,outs,nproc,iproc,infocode)
!      ncount_bigdft=ncount_bigdft+1
!      call vcopy(3 * outs%fdim, outs%fxyz(1,1), 1, fxyz(1,1,nhist), 1)
!      etotp=outs%energy
      detot=etotp-etotold

      if(debug.and.iproc==0)write(100,'(a,i6,2(1x,e21.14),1x,5(1x,e10.3),xi0)')&
            'SBFGS it,etot,etotold,Detot,fnrm,fnrmp/fnrm,dnrm/fnrm,beta,ndim',&
             it-1,etotp,etotold,detot,fnrm,sqrt(ts)/fnrm,sqrt(tt)/fnrm,beta,ndim


      call fnrmandforcemax(fxyz(1,1,nhist),fnrm,fmax,nat)
      fnrm=sqrt(fnrm)

      if (iproc == 0) then
         write(fn4,'(i4.4)') ncount_bigdft
         write(comment,'(a,1pe10.3)')'SBFGS:fnrm= ',fnrm
         call write_atomic_file(trim(runObj%inputs%dir_output)//'posout_'//fn4,&
              outs%energy,runObj%atoms%astruct%rxyz,runObj%atoms%astruct%ixyz_int,&
              runObj%atoms,trim(comment),forces=outs%fxyz)
      endif

      if (fmax < 3.d-1) call updatefluctsum(outs%fnoise,fluct)
      cosangle=-dot_double(3*nat,fxyz(1,1,nhist),1,dd(1,1),1)/&
              sqrt(dot_double(3*nat,fxyz(1,1,nhist),1,fxyz(1,1,nhist),1)*&
              dot_double(3*nat,dd(1,1),1,dd(1,1),1))


   
      if (detot.gt.maxrise .and. beta > 1.d-1*betax) then !
         if (debug.and.iproc==0) write(100,'(a,i0,1x,e9.2)') "WARN: it,detot", it,detot
         if (debug.and.iproc==0) write(16,'(a,i0,4(xe9.2))') &
             "WARNING GEOPT_SBFGS: Prevent energy to rise by more than maxrise: it,maxrise,detot,beta,1.d-1*betax ",&
             it,maxrise,detot,beta,1.d-1*betax
         if (iproc==0.and.verbosity > 0) then
            !avoid space for leading sign (numbers are positive, anyway)
            write(cdmy8,'(es8.1)')abs(maxd)
            write(cdmy9_1,'(es9.2)')abs(displr)
            write(cdmy9_2,'(es9.2)')abs(displp)
            write(cdmy9_3,'(es9.2)')abs(beta)
   
            write(16,'(i5,1x,i5,2x,a10,2x,1es21.14,2x,es9.2,es11.3,3es10.2,2x,a6,a8,xa4,i3.3,xa5,a7,2(xa6,a8))') &
             ncount_bigdft,it,'GEOPT_SBFGS',etotp,detot,fmax,fnrm,fluct*runObj%inputs%frac_fluct,fluct, &
             'beta=',trim(adjustl(cdmy9_3)),'dim=',ndim,'maxd=',trim(adjustl(cdmy8)),'dsplr=',trim(adjustl(cdmy9_1)),'dsplp=',trim(adjustl(cdmy9_2))
            call yaml_open_map('Geometry')
               call yaml_map('Ncount_BigDFT',ncount_bigdft)
               call yaml_map('Geometry step',it)
               call yaml_map('Geometry Method','GEOPT_SBFGS')
               call yaml_map('ndim',ndim)
               call yaml_map('etot', etotp,fmt='(1pe21.14)')
               call yaml_map('detot',detot,fmt='(1pe21.14)')
               call yaml_map('fmax',fmax,fmt='(1pe21.14)')
               call yaml_map('fnrm',fnrm,fmt='(1pe21.14)')
               call yaml_map('beta',beta,fmt='(1pe21.14)')
               call geometry_output(fmax,fnrm,fluct)
            call yaml_close_map()
         end if
    
         if(ncount_bigdft >= nit)then!no convergence within ncount_cluster_x energy evaluations
            !following copy of rxyz(1,1,nhist-1) to runObj is necessary for returning to the caller
            !the energies and coordinates used/obtained from/in the last ACCEPTED iteration step
            !(otherwise coordinates of last call to call_bigdft would be returned)
            call vcopy(3 * runObj%atoms%astruct%nat, rxyz(1,1,nhist-1), 1,runObj%atoms%astruct%rxyz(1,1), 1)
            goto 900  !sbfgs will return to caller the energies and coordinates used/obtained from the last ACCEPTED iteration step
         endif

         !beta=min(.50_gp*beta,betax)
         !beta=max(.50_gp*beta,1.d-1*betax)
         beta=.50_gp*beta
         if (debug.and.iproc==0) write(100,'(a,1x,e9.2)') 'WARNING GEOPT_SBFGS: beta reset ',beta
         ndim=0
         if(.not.steep)then
            do iat=1,nat
               rxyz(1,iat,0)=rxyz(1,iat,nhist-1)
               rxyz(2,iat,0)=rxyz(2,iat,nhist-1)
               rxyz(3,iat,0)=rxyz(3,iat,nhist-1)
   
               fxyz(1,iat,0)=fxyz(1,iat,nhist-1)
               fxyz(2,iat,0)=fxyz(2,iat,nhist-1)
               fxyz(3,iat,0)=fxyz(3,iat,nhist-1)
            enddo
            nhist=1
         endif
         goto  500
      endif

      displp=displp+tt
      if (iproc==0.and.verbosity > 0) then
         !avoid space for leading sign (numbers are positive, anyway)
         write(cdmy8,'(es8.1)')abs(maxd)
         write(cdmy9_1,'(es9.2)')abs(displr)
         write(cdmy9_2,'(es9.2)')abs(displp)
         write(cdmy9_3,'(es9.2)')abs(beta)

         write(16,'(i5,1x,i5,2x,a10,2x,1es21.14,2x,es9.2,es11.3,3es10.2,2x,a6,a8,xa4,i3.3,xa5,a7,2(xa6,a8))') &
          ncount_bigdft,it,'GEOPT_SBFGS',etotp,detot,fmax,fnrm,fluct*runObj%inputs%frac_fluct,fluct, &
          'beta=',trim(adjustl(cdmy9_3)),'dim=',ndim,'maxd=',trim(adjustl(cdmy8)),'dsplr=',trim(adjustl(cdmy9_1)),'dsplp=',trim(adjustl(cdmy9_2))
         call yaml_open_map('Geometry')
            call yaml_map('Ncount_BigDFT',ncount_bigdft)
            call yaml_map('Geometry step',it)
            call yaml_map('Geometry Method','GEOPT_SBFGS')
            call yaml_map('ndim',ndim)
            call yaml_map('etot', etotp,fmt='(1pe21.14)')
            call yaml_map('detot',detot,fmt='(1pe21.14)')
            call yaml_map('fmax',fmax,fmt='(1pe21.14)')
            call yaml_map('fnrm',fnrm,fmt='(1pe21.14)')
            call yaml_map('beta',beta,fmt='(1pe21.14)')
            call geometry_output(fmax,fnrm,fluct)
         call yaml_close_map()
      end if

      etot    = etotp
      etotold = etot
      !copy outs_datatype
!      call copy_global_output(outs,outsIO)

      if(detot .gt. maxrise)then
         if (iproc==0) write(16,'(a,i0,4(xe9.2))') &
             "WARNING GEOPT_SBFGS: Allowed energy to rise by more than maxrise: it,maxrise,detot,beta,1.d-1*betax ",&
             it,maxrise,detot,beta,1.d-1*betax
      endif


!      if (fnrm.le.fnrmtol) goto 1000
      call convcheck(fmax,fluct*runObj%inputs%frac_fluct,runObj%inputs%forcemax,icheck)
      if(icheck>5)then
         goto 1000
      endif

      if(ncount_bigdft >= nit)then!no convergence within ncount_cluster_x energy evaluations
            goto 900  !sbfgs will return to caller the energies and coordinates used/obtained from the last accepted iteration step
      endif
   
      if (cosangle.gt..200_gp) then
         beta=beta*1.100_gp
      else
         beta=max(beta*.850_gp,betax)
      endif
   
      if (debug.and.iproc==0) write(100,*) 'cosangle ',cosangle,beta

      ! calculate norms
      do i=1,nhist
           rnorm(i)=0.d0
              do iat=1,nat
                 do l=1,3
                    rnorm(i)=rnorm(i) + (rxyz(l,iat,i)-rxyz(l,iat,i-1))**2
                 enddo
              enddo
     rnorm(i)=1.d0/sqrt(rnorm(i))
     !rnorm(i)=1.d0*sqrt(dble(i))/sqrt(rnorm(i))
      !rnorm(i) =
      !(0.5d0*(1.d0+(((dble(i)-dble(nhist))/(dble(nhist)*0.25d0))/sqrt(1.d0+((dble(i)-dble(nhist))/(dble(nhist)*0.25d0))**2))))/sqrt(rnorm(i))
      !rnorm(i) = (0.5d0*(1.d0+(((dble(i)-0.2d0*dble(nhist))/(dble(nhist)*0.25d0))/sqrt(1.d0+((dble(i)-0.2d0*dble(nhist))/(dble(nhist)*0.25d0))**2))))/sqrt(rnorm(i))

      enddo

   
      !find linear dependencies via diagonalization of overlap matrix   
      !build overlap matrix:
      do i=1,nhist
         do j=1,nhist
            aa(i,j)=0.0_gp
            do iat=1,nat
               do l=1,3
                  aa(i,j)=aa(i,j)&
                  +(rxyz(l,iat,i)-rxyz(l,iat,i-1))&
                  *(rxyz(l,iat,j)-rxyz(l,iat,j-1))
               enddo
            enddo
            aa(i,j)=aa(i,j)*rnorm(i)*rnorm(j)
         enddo
      enddo
   
      !diagonalize overlap matrix:
      call dsyev('V',"L",nhist,aa,nhistx,eval,work,lwork,info)
      if (info.ne.0) then 
         if (debug.and.iproc==0) write(100,*) ' Over ', info
         stop 'info geo after first dsyev aa'
      endif
      if(debug.and.iproc==0)then
         do i=1,nhist
            write(100,*) "Overl ",i,eval(i)
         enddo
      endif
      
      do idim=1,nhist
         do iat=1,nat
            do l=1,3
               rr(l,iat,idim)=0.0_gp
               ff(l,iat,idim)=0.0_gp
            enddo
         enddo
      enddo
     
      !generate significant orthogonal subspace 
      ndim=0
      do idim=1,nhist
         !remove linear dependencies by using the overlap-matrix eigenvalues:
         if (eval(idim)/eval(nhist).gt.cutoffRatio) then    ! HERE
            ndim=ndim+1
      
            do jdim=1,nhist
               do iat=1,nat
                  do l=1,3
                     rr(l,iat,ndim)=rr(l,iat,ndim)+aa(jdim,idim)*rnorm(jdim)&
                                 *(rxyz(l,iat,jdim)-rxyz(l,iat,jdim-1))
                     ff(l,iat,ndim)=ff(l,iat,ndim)-aa(jdim,idim)*rnorm(jdim)&
                                 *(fxyz(l,iat,jdim)-fxyz(l,iat,jdim-1))
                  enddo
               enddo
            enddo
      
            do iat=1,nat
               do l=1,3
                  rr(l,iat,ndim)=rr(l,iat,ndim)/sqrt(abs(eval(idim)))
                  ff(l,iat,ndim)=ff(l,iat,ndim)/sqrt(abs(eval(idim)))
               enddo
            enddo
         endif
      enddo
      if (debug.and.iproc==0) write(100,'(a,i3)') "ndim= ",ndim
      
      ! Hessian matrix in significant orthogonal subspace
      do i=1,ndim
         do j=1,ndim
            aa(i,j)=0.0_gp
            do iat=1,nat
               do l=1,3
                  aa(i,j)=aa(i,j) + .50_gp*(rr(l,iat,i)*ff(l,iat,j)&
                                  +rr(l,iat,j)*ff(l,iat,i))
               enddo
            enddo
         enddo
      enddo
      
      call dsyev('V',"L",ndim,aa,nhistx,eval,work,lwork,info)
      if (info.ne.0) then 
         write(*,*) 'ERROR: info after 2n dsyev aa in sbfgs', info
         stop 'info after 2nd dsyev aa in sbfgs'
      endif
      
      ! calculate eigenvectors
      do i=1,ndim
         do iat=1,nat
            do l=1,3
               rrr(l,iat,i)=0.0_gp
               fff(l,iat,i)=0.0_gp
            enddo
         enddo
      enddo
      
      do i=1,ndim
         tt=0.0_gp
         do j=1,ndim
            do iat=1,nat
               do l=1,3
                  rrr(l,iat,i)=rrr(l,iat,i) + aa(j,i)*rr(l,iat,j)
                  fff(l,iat,i)=fff(l,iat,i) + aa(j,i)*ff(l,iat,j)
               enddo
            enddo
         enddo
         do iat=1,nat
            do l=1,3
               tt=tt+(fff(l,iat,i)-eval(i)*rrr(l,iat,i))**2
            enddo
         enddo
         !residuue according to Weinstein criterion
         res(i)=sqrt(tt)
         if (debug.and.iproc==0) write(100,'(a,i3,e14.7,1x,e12.5,2(1x,e9.2))') 'EVAL,RES '&
                                                             ,i,eval(i),res(i)
      enddo
   enddo!end main loop

900 continue

   !if code gets here, it failed
   if(debug.and.iproc==0) write(100,*) it,etot,fnrm
   if(iproc==0) write(16,'(a,3(xi0))') &
       "WARNING GEOPT_SBFGS: SBFGS not converged: it,ncount_bigdft,ncount_cluster_x: ", &
       it,ncount_bigdft,runObj%inputs%ncount_cluster_x
!   stop "No convergence "
   fail=.true.
   goto 2000

1000 continue!converged successfully
   
   if(iproc==0) write(16,'(2(a,xi0))') "SBFGS converged at iteration ",it,". Needed bigdft calls: ",ncount_bigdft
   if(iproc==0)  call yaml_map('Iterations when SBFGS converged',it)
   fail=.false.
   
!   etot=etotp
!   do iat=1,nat
!      do l=1,3
!         wpos(l,iat)= rxyz(l,iat,nhist)
!         fout(l,iat)= fxyz(l,iat,nhist)
!      enddo
!   enddo
2000 continue
!deallocations
   call f_free(rxyz)
   call f_free(fxyz)
   call f_free(aa)
   call f_free(eval)
   call f_free(res)
   call f_free(rnorm)
   call f_free(ff)
   call f_free(rr)
   call f_free(dd)
   call f_free(fff)
   call f_free(rrr)
   call f_free(scpr)
!   call deallocate_global_output(outs)
end subroutine
end module
