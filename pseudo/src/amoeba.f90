!> @file
!! atomic program for generating and optimizing HGH pseudo-potentials.
!! @author
!!    Alex Willand, under the supervision of Stefan Goedecker
!!    gpu accelerated routines by Raffael Widmer
!!    parts of this program were based on the fitting program by matthias krack
!!    http://cvs.berlios.de/cgi-bin/viewcvs.cgi/cp2k/potentials/goedecker/pseudo/v2.2/


!> Minimizes the penalty function (subroutine) using
!! a simplex downhill method (amoeba)
!! @note cwh: if amoeba doesn't improve within ntrymax iterations amoeba also returns
subroutine amoeba(p,y,ndim,ftol,iter,itmax,namoeb, iproc,nproc, ntrymax, energ,verbose)
   
   use pseudovars
   use penaltyvars
   
   implicit none
   !Arguments
   integer, intent(in) :: ndim, iproc,nproc,itmax,ntrymax
   integer, intent(inout) :: namoeb
   real(kind=8), intent(in) :: ftol,energ
   logical, intent(out) :: verbose
   real(kind=8), dimension(ndim,ndim+1), intent(out) :: p
   real(kind=8), dimension(ndim+1), intent(out) :: y
   !Local variables
   real(kind=8), parameter :: alpha=1.0d0,beta=0.5d0,gamma=2.0d0
   logical :: lexit,lnext
   real(kind=8), dimension(ndim) :: pr,prr,pbar
   real(kind=8) :: rtol,ypr,yprr
   integer :: mpts,iter,ntrycount,ilo,ihi,inhi,i,ierr,iloold,j
   include 'mpif.h'
   
   ! print*,'entered amoeba with ',ndim,itmax,ftol,ntrymax
   mpts=ndim+1
   iter=0
   ntrycount = 0
   iloold = 1
   
   ! ================== here begins the loop for the simplex ======================
   loop_simplex: do

      !first, get the currently highest, second highest and lowest vertices of the simplex
      ilo=1
      if (y(1).gt.y(2))then
         ihi=1
         inhi=2
      else
         ihi=2
         inhi=1
      endif
      do i=1,mpts
         if (y(i).lt.y(ilo)) ilo=i
         if (y(i).gt.y(ihi))then
            inhi=ihi
            ihi=i
         else if (y(i).gt.y(inhi))then
            if (i.ne.ihi) inhi=i
         endif
         
      end do
      
      ! then do various checks for exit conditions
      
      ! cwh
      if (ilo .ne. iloold) ntrycount = 0
      ! write(6,*)'debug: ilo, ntrycount,ntrymax',ilo, ntrycount,ntrymax
      iloold = ilo
      
      ! rtol=min(y(ihi),(y(ihi)-y(ilo))/y(ilo)**2)
      rtol=min(y(ihi),(y(ihi)-y(ilo))/y(ilo))
      if (rtol < ftol) then
         ! check
         write(6,'(1x,a,1pe20.10,a,1pe20.10)') 'amoeba converged:',rtol,' <',ftol
         ! write(6,*) 'values at the edges of the simplex:'
         ! write(6,'(40e15.8)') y
         do i=1,ndim
            if (y(i).lt.y(ilo)) write(6,*) 'warning ilo not lowest'
         enddo
      endif
      ! if (mod(iter,100).eq.0) write(6,*) 'iter=',iter
      
      if (iter.gt.itmax) then
         ! write(6,*) 'warning: no convergence in amoeba'
         write(6,'(1x,a)') '--- simplex done, max no of iterations reached.'
         ! ftol=10.d0*ftol
         ! write(6,*) 'ftol set to ',ftol
         ! write(6,*) 'values at the edges of the simplex:'
         ! write(6,'(40e15.8)') y
         ! write(6,*)'ilo:',ilo
         do i=1,ndim
            if (y(i).lt.y(ilo)) write(6,'(1x,a)') 'warning ilo not lowest'
         enddo
      end if
      
      if (ntrycount.ge.ntrymax) then
         write(6,*) 'no improvement during',ntrycount, 'iterations, tolerance is exceeded.'
         ! write(6,*) 'warning: no improvement in amoeba for the last', ntrycount,'steps'
      endif
      
      ! every 10th step, call a barrier and check for exit/next files
      if (mod(iter,10).eq.0) then
         ! a barrier for synchronisation cannot be wrong here
         if (nproc > 1) call mpi_barrier(mpi_comm_world,ierr)  
         inquire ( file = 'next', exist = lnext )
         inquire ( file = 'exit', exist = lexit )
         if (lnext) then 
            if (nproc > 1) call mpi_barrier(mpi_comm_world,ierr)  
            write(6,*) 'the file next is present, aborting this fit' 
            if (iproc == 0)then
               open(unit=99,file='next')
               close(unit=99, status='delete')
            end if
         endif
         if (lexit) then
            ! set namoeb to something smaller than iiter so that the main programs cycle loop exits 
            namoeb=-namoeb
            if (nproc > 1) call mpi_barrier(mpi_comm_world,ierr)  
            write(6,*) 'the file exit is present, aborting fits' 
            if (iproc.eq.0)then
               open(unit=99,file='exit')
               close(unit=99, status='delete')
            end if
         endif
      endif
      
      
      ! several exit conditions to leave amoeba, but always the same actions:
      !         pack the lowest vertex in the current psppar
      !         call penalty to make sure the latest call to gatom was indeed
      !         for the lowest vertex and therefore the packed psppar
      !         call penalty rather than gatom to give some more information
      !         about the excitation energies and the softness when leaving
      !         the simplex cycle
      
      if ( (iter.gt.itmax) .or. lexit .or. lnext .or. (ntrycount.ge.ntrymax) .or. (rtol.lt.ftol) ) then
         
         ! another call to penalty with the verbose flag yields further information
         write(6,'(1x,a)')     '______________'
         write(6,'(1x,a,/,/)') 'leaving amoeba'
         verbose=.true.
         ! write(6,*)'verbose=',verbose,' inquire penalty information'
         
         call penalty(energ,verbose,p(1,ilo),y(ilo), iproc,nproc,iter,y(ilo))
         
         ! unpack variables from p(ilo) into psppar variables
         verbose=.false.
         call ppack (verbose, p(1,ilo), "unpack") 
         ! verbose=.false.
         return
      endif
      
      ! cwh
      ntrycount = ntrycount +1
      if (mod(ntrycount,max(1,ntrymax/5)) == 0) &
         write(6,*)'no improvement during',ntrycount, 'iterations, tolerance is',ntrymax
      
      ! ================== here are the actual simplex downhill moves ======================
      iter=iter+1
      ! get the simplex centre
      do j=1,ndim
         pbar(j)=0.d0
      end do
      do i=1,mpts
         if (i.ne.ihi)then
            do j=1,ndim
               pbar(j)=pbar(j)+p(j,i)
            end do
         endif
      end do
      
      ! reflect the highest vertex at the centre -> pr
      do j=1,ndim
         pbar(j)=pbar(j)/ndim
         pr(j)=(1.d0+alpha)*pbar(j)-alpha*p(j,ihi)
      end do
      call penalty(energ,verbose,pr,ypr,  &
           iproc,nproc,iter,y(ilo))
      if (ypr.le.y(ilo))then
         ! new lowest, can we go further in that direction?
         
         ! test a bigger step for the reflecting 
         do  j=1,ndim
            prr(j)=gamma*pr(j)+(1.d0-gamma)*pbar(j)
         end do
         call penalty(energ,verbose,prr,yprr,  &
              iproc,nproc,iter,y(ilo))
         if (yprr.lt.y(ilo))then
            ! this got even better, so accept prr, discarding the highest vertex 
            do  j=1,ndim
               p(j,ihi)=prr(j)
            end do
            ! if (mod(iter,10).eq.0) write(6,'(a,i5,a,e15.7,a,i2,e15.7)') &
            !    'iter',iter,' found',yprr,' reject',ihi,y(ihi)
            y(ihi)=yprr
         else
            ! pr is the best we have, so accept it, discarding the highest vertex 
            do j=1,ndim
               p(j,ihi)=pr(j)
            end do
            ! if (mod(iter,10).eq.0) write(6,'(a,i5,a,e15.7,a,i2,e15.7)') &
            !    'iter',iter,' found',ypr,' reject',ihi,y(ihi)
            y(ihi)=ypr
         endif
      else if (ypr.ge.y(inhi))then
         ! the reflected vertex is not lower than the second highest vertex
         if (ypr.lt.y(ihi))then
            ! if reflecting improved the highest vertex, update it
            do j=1,ndim
               p(j,ihi)=pr(j)
            end do
            ! if (mod(iter,10).eq.0) write(6,'(a,i5,a,e15.7,2(a,i2,e15.7))') &
            !     'iter',iter,' found',ypr,' reject',ihi,y(ihi), ' best:',ilo,y(ilo)
            y(ihi)=ypr
         endif
         ! try to contract the highest vertex to the centre
         do j=1,ndim
            prr(j)=beta*p(j,ihi)+(1.d0-beta)*pbar(j)
         end do
         call penalty(energ,verbose,prr,yprr,  &
              iproc,nproc,iter,y(ilo))
         if (yprr.lt.y(ihi))then
            ! if contracting improved the highest vertex, save it
            do j=1,ndim
               p(j,ihi)=prr(j)
            end do
            ! if (mod(iter,10).eq.0) write(6,'(a,i5,a,e15.7,2(a,i2,e15.7))') &
            !    'iter',iter,' found',yprr,' reject',ihi,y(ihi),' best:',ilo,y(ilo)
            y(ihi)=yprr
         else
            ! all moves of the highest vertex failed to improof the penalty,
            ! so we shrink the entire simplex towards the lowest vertex
            do i=1,mpts
               if (i.ne.ilo)then
                  do j=1,ndim
                     pr(j)=0.5d0*(p(j,i)+p(j,ilo))
                     p(j,i)=pr(j)
                  end do
                  call penalty(energ,verbose,pr,ypr,  &
                       iproc,nproc,iter,y(ilo))
               endif
            end do
         endif
      else
         do j=1,ndim
            p(j,ihi)=pr(j)
         end do
         ! if (mod(iter,10).eq.0) write(6,'(a,i5,a,e15.7,2(a,i2,e15.7))') &
         !    'iter',iter,' found',ypr,' reject',ihi,y(ihi), ' best:',ilo,y(ilo)
         y(ihi)=ypr
      endif

   end do loop_simplex !< Final loop over vertices of the simplex

end subroutine amoeba
