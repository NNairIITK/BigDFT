!> @file
!! Wannier constructor
!! @author
!!    Copyright (C) 2011-2012 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


program WaCo

   use module_base
   use module_types
   use module_interfaces, except_this_one => writeonewave
   use Poisson_Solver
   implicit none
   character :: filetype*4,outputype*4
   type(locreg_descriptors) :: Glr
   type(orbitals_data) :: orbs,orbsw,orbsv
   type(atoms_data) :: atoms
   type(input_variables) :: input
   type(workarr_sumrho) :: w
   type(communications_arrays), target :: commsw
   real(gp), parameter :: b2a=0.5291772108_dp
   real :: tcpu0,tcpu1
   real(gp) :: tel
   real(gp), dimension(3) :: shift,CM
   real(gp) :: dist,rad,sprdfact,sprddiff,enediff
   integer :: iproc, nproc, nproctiming, i_stat, nelec, ierr, npsidim
   integer :: nvirtu,nvirtd,nrpts
   integer :: NeglectPoint, CNeglectPoint
   integer :: ncount0,ncount1,ncount_rate,ncount_max,iat,iformat
   integer :: iwann, iiwann, iband, nwann, nband, plotwann
   integer :: iw1, iw2, ifile
   integer :: nsprd,ndiag, nwannCon
   character(len=*), parameter :: subname='WaCo'
   character(len=4) :: num, units
   integer, allocatable :: wann_list(:), Zatoms(:,:), ncenters(:), types(:,:)
   real(gp), dimension(:,:), pointer :: rxyz, rxyz_old, cxyz,rxyz_wann
   real(gp), dimension(:,:), allocatable :: radii_cf
   real(gp), allocatable :: sprd(:), eigen(:,:), proj(:,:), projC(:,:),distw(:),charge(:),prodw(:),wannocc(:)
   real(wp), allocatable :: psi(:,:),wann(:),wannr(:)
   real(wp), allocatable :: ham(:,:,:),hamr(:,:,:)
   real(wp), allocatable :: diag(:,:),diagT(:)
   integer, dimension(:), pointer :: buf
   character(len=60) :: radical, filename
   character(len=8) :: string
   logical :: notocc, bondAna,Stereo,hamilAna,WannCon
   integer, dimension(:), allocatable :: ConstList
   integer, allocatable :: nfacets(:),facets(:,:,:),vertex(:,:,:), l(:), mr(:)
   real(gp), dimension(3) :: refpos, normal
   real(kind=8), allocatable :: umn(:,:), rho(:,:), rhoprime(:,:),amn(:,:),tmatrix(:,:)
   integer :: i, j, k, i_all
   character(len=16) :: seedname
   integer :: n_occ, n_virt, n_virt_tot, nproj,nband_old,nkpt_old 
   logical :: w_unk, w_sph, w_ang, w_rad, pre_check
   integer, allocatable, dimension (:) :: virt_list
   logical :: file_exist,idemp
!   real(gp),allocatable :: cxyz2(:,:) !debug only
!   integer, allocatable :: list(:)    !debug only
   interface
      subroutine read_inter_header(iproc,seedname, filetype, n_occ, pre_check, n_virt_tot, n_virt, w_unk, w_sph, w_ang, w_rad)
        implicit none
        integer, intent(in) :: iproc
        character, intent(out) :: seedname*16, filetype*4
        integer, intent(out) :: n_occ, n_virt, n_virt_tot
        logical, intent(out) :: w_unk, w_sph, w_ang, w_rad, pre_check
      end subroutine read_inter_header
      subroutine read_inter_list(iproc,n_virt, virt_list)
        implicit none
        integer, intent(in) :: n_virt,iproc
        integer, dimension(n_virt), intent(out) :: virt_list
      end subroutine read_inter_list
      subroutine scalar_kmeans_diffIG(iproc,nIG,crit,nel,vect,string,nbuf,buf)
        implicit none
        integer, intent(in) :: nel,nIG,iproc
        real(kind=8),intent(in) :: crit
        real(kind=8), dimension(nel),intent(in) :: vect
        character(len=*),intent(in) :: string
        integer, intent(out) :: nbuf
        integer, dimension(:), pointer :: buf
      end subroutine scalar_kmeans_diffIG
   end interface

   ! Start MPI in parallel version
   !in the case of MPIfake libraries the number of processors is automatically adjusted
   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

   call memocc_set_memory_limit(memorylimit)

   ! Read a possible radical format argument.
   call get_command_argument(1, value = radical, status = i_stat)
   if (i_stat > 0) then
      write(radical, "(A)") "input"
   end if

   if (input%verbosity > 2) then
      nproctiming=-nproc !timing in debug mode                                                                                                                                                                 
   else
      nproctiming=nproc
   end if

   call timing(nproctiming,'WaCo_time.prc','IN')

   call cpu_time(tcpu0)
   call system_clock(ncount0,ncount_rate,ncount_max) 
   call timing(iproc,'Precondition  ','ON')   

   !###################################################################
   ! Read input files and initialise the variables for the wavefunctions
   !###################################################################
   call standard_inputfile_names(input,radical)

   call Waco_input_variables(iproc,trim(radical)//'.waco',nband,nwann,bondAna,Stereo,hamilAna,WannCon,&
        outputype,nwannCon,refpos,units,sprdfact,sprddiff,enediff)

   allocate(ConstList(nwannCon))
   call read_input_waco(trim(radical)//'.waco',nwannCon,ConstList) 

   call read_input_variables(iproc,'posinp',input, atoms, rxyz)

   if (iproc == 0) call print_general_parameters(nproc,input,atoms)

   allocate(radii_cf(atoms%ntypes,3+ndebug),stat=i_stat)
   call memocc(i_stat,radii_cf,'radii_cf',subname)

   call system_properties(iproc,nproc,input,atoms,orbs,radii_cf,nelec)

   ! Determine size alat of overall simulation cell and shift atom positions
   ! then calculate the size in units of the grid space
   call system_size(iproc,atoms,rxyz,radii_cf,input%crmult,input%frmult,input%hx,input%hy,input%hz,&
        Glr,shift)

   ! Create wavefunctions descriptors and allocate them inside the global locreg desc.
   call createWavefunctionsDescriptors(iproc,input%hx,input%hy,input%hz,&
        atoms,rxyz,radii_cf,input%crmult,input%frmult,Glr)

   ! don't need radii_cf anymore
   i_all = -product(shape(radii_cf))*kind(radii_cf)
   deallocate(radii_cf,stat=i_stat)
   call memocc(i_stat,i_all,'radii_cf',subname)

   !#################################################################
   ! Read the input.inter file
   !#################################################################
   call read_inter_header(iproc,seedname, filetype, n_occ, pre_check, n_virt_tot, n_virt, w_unk, w_sph, w_ang, w_rad)
   allocate(virt_list(n_virt),stat=i_stat)
   call memocc(i_stat,virt_list,'virt_list',subname)
   if (n_virt .ne. 0) then
      call read_inter_list(iproc, n_virt, virt_list)
   end if 

   !#####################################################################
   ! Read the Umn file
   !#####################################################################
!   nwann = n_occ + n_virt             !for now, only square matrices
!   nband = n_occ + n_virt

   if ((nband .ne. n_occ+n_virt) .and. iproc == 0) then
      write(*,*) 'Number of bands in the .waco file : ',nband
      write(*,*) 'not equal to the number of bands used in .inter file:', n_occ+n_virt
      call mpi_finalize(ierr)
      stop
   end if

   allocate(umn(nwann,nband),stat=i_stat)
   call memocc(i_stat,umn,'umn',subname)
   call read_umn(iproc,nwann,nband,seedname,umn)
   call timing(iproc,'Precondition  ','OF')

   !###################################################################
   ! Constructing the density matrix
   !###################################################################

   if(iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '!   Constructing density matrix    !'
      write(*,*) '!==================================!'
      write(*,*)
   end if

   allocate(rho(nwann,nwann),stat=i_stat)
   call memocc(i_stat,rho,'rho',subname)
   call to_zero(nwann*nwann, rho(1,1))

   do i=1,nwann
      do j=1,nwann
         do iband=1,nband
            rho(i,j) = rho(i,j) + umn(i,iband)*umn(j,iband)
          end do
      end do
   end do

!!   if(iproc==0) then
!!      do i=1,nwann
!!         do j=1,nwann
!!            write(*,*) i , j ,rho(i,j)
!!         end do
!!      end do
!!   end if

   if(iproc==0) then
      write(*,*) '!===================================!'
      write(*,*) '!Constructing density matrix : DONE !'
      write(*,*) '!===================================!'
      write(*,*)
   end if

   !##################################################################
   ! Check idempotence of density matrix
   !##################################################################
   allocate(rhoprime(nwann,nwann),stat=i_stat)
   call memocc(i_stat,rhoprime,'rhoprime',subname)
   call to_zero(nwann*nwann, rhoprime(1,1))
   
   idemp = .true.
   do i = 1, nwann
      do j = 1, nwann
         do k = 1, nwann
            rhoprime(i,j) = rhoprime(i,j) + rho(i,k)*rho(k,j)
         end do
         if(abs(rhoprime(i,j) - rho(i,j))>1.0d-5) then
           write(*,*) 'Not indempotent',i,j,rhoprime(i,j), rho(i,j)
           idemp = .false.
         end if
      end do
   end do

   if(idemp .eqv. .false.) then
     stop 'Density matrix not idempotent'
   end if

   i_all = -product(shape(rhoprime))*kind(rhoprime)
   deallocate(rhoprime,stat=i_stat)
   call memocc(i_stat,i_all,'rhoprime',subname)
   i_all = -product(shape(rho))*kind(rho)
   deallocate(rho,stat=i_stat)
   call memocc(i_stat,i_all,'rho',subname)

   !#########################################################
   ! Read seedname.dat to make the Wannier plot list
   !#########################################################
   if(iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '!     Reading .dat :               !'
      write(*,*) '!==================================!'
   end if

   inquire(file=trim(seedname)//'.dat',exist=file_exist)
   if (.not. file_exist) then
      if (iproc==0) then
         write(*,'(A,1x,A)') 'ERROR : Input file,',trim(seedname)//'.dat, not found !'
         write(*,'(A)') 'CORRECTION: Create or give correct input file.'
      end if
      call mpi_finalize(i)
      stop
   end if

   !wann_list will contain the list of occupied Wannier functions
   allocate(wann_list(nwann),stat=i_stat)
   call memocc(i_stat,wann_list,'wann_list',subname)
   wann_list = 0

   open(11, file=trim(seedname)//'.dat', status='OLD')

   ! Check if the line is commented
   plotwann = 0
   do i = 1 ,nwann
      read(11,*) string
      if(VERIFY(string, ' #', .false.) == 0) cycle
      plotwann = plotwann +1
      wann_list(plotwann) = i
   end do

   close(11)

   if(iproc == 0) then
      write(*,'(A)') 'Occupied/Plotted Wannier functions are:'
      do i=1, plotwann, 6
         write(*,'(6(2x,I4))') (wann_list(j), j=i,min(i+5,size(wann_list)))
      end do
   end if

   if(iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '!     Reading .dat : DONE          !'
      write(*,*) '!==================================!'
   end if

   !########################################################
   ! Bonding analysis
   !########################################################
  if(bondAna) then
     if(iproc==0) then
         write(*,*) '!==================================!'
         write(*,*) '!     Bonding Analysis :           !'
         write(*,*) '!==================================!'
     end if


      !wann_list will contain the list of occupied Wannier functions
      allocate(cxyz(3,plotwann),stat=i_stat)
      call memocc(i_stat,cxyz,'cxyz',subname)
      allocate(rxyz_wann(3,atoms%nat),stat=i_stat)
      call memocc(i_stat,rxyz_wann,'rxyz_wann',subname)
      allocate(sprd(plotwann+1),stat=i_stat)
      call memocc(i_stat,sprd,'sprd',subname)
      allocate(Zatoms(atoms%nat,plotwann),stat=i_stat)
      call memocc(i_stat,Zatoms,'Zatoms',subname)  
      allocate(ncenters(plotwann),stat=i_stat)
      call memocc(i_stat,ncenters,'ncenters',subname)

      call read_centers(iproc,nwann,plotwann,atoms%nat,seedname,wann_list,cxyz,rxyz_wann,.true.,sprd)

      ! Now calculate the bonding distances and ncenters
      if(iproc == 0) write(*,'(2x,A)') 'Number of atoms associated to the WFs'
      if(iproc == 0) write(*,'(3x,A,x,A,4x,A,3x,A,3x,A)') 'WF','OWF','Spr(ang^2)','Nc','Atom numbers:'
      Zatoms = 0
      do iwann = 1, plotwann
         ncenters(iwann) = 0
         iat = 0
         do i = 1, atoms%nat
            dist = (rxyz_wann(1,i)-cxyz(1,iwann))**2 + (rxyz_wann(2,i)-cxyz(2,iwann))**2 + (rxyz_wann(3,i)-cxyz(3,iwann))**2
            if (dist <= sprdfact * sprd(iwann)) then    !for normal distribution: 1=68%, 1.64=80%, 3=94%
               ncenters(iwann) = ncenters(iwann) +1
               iat = iat +1
               Zatoms(iat,iwann) = i
            end if
         end do
         if(iproc == 0) then
           write(*,'(2I4,F14.6,2x,I4,6(2x,I4))') wann_list(iwann),iwann, sprd(iwann), ncenters(iwann),&
                (Zatoms(i_all,iwann),i_all=1,iat)
         end if
      end do

     ! Calculate occupation of the wannier functions
     allocate(wannocc(nwann),stat=i_stat)
     call memocc(i_stat,wannocc,'wannocc',subname)
     wannocc = 0.0_dp
     do iwann = 1, nwann
        do iband = 1, n_occ
            wannocc(iwann) = wannocc(iwann) + umn(iwann,iband)**2
        end do
     end do
     print *,'total number of electrons: ',2*sum(wannocc)

      allocate(distw(maxval(ncenters)),stat=i_stat)
      call memocc(i_stat,distw,'distw',subname)
      allocate(charge(atoms%nat),stat=i_stat)
      call memocc(i_stat,charge,'charge',subname)
      allocate(prodw(maxval(ncenters)),stat=i_stat)
      call memocc(i_stat,prodw,'prodw',subname)

      ! Calculate "Wannier charge" = 2 *occupation* product(d_all_atoms__associated_with_this_wannier_except_this_one) / d_total
      charge = 0.0_dp
      do iwann = 1, plotwann
         distw = 0.0_dp
         do iat = 1, ncenters(iwann)
            distw(iat) = sqrt((rxyz_wann(1,Zatoms(iat,iwann))-cxyz(1,iwann))**2 +&
                        (rxyz_wann(2,Zatoms(iat,iwann))-cxyz(2,iwann))**2 +&
                        (rxyz_wann(3,Zatoms(iat,iwann))-cxyz(3,iwann))**2)
         end do
         prodw = 0.0_dp
         do iat = 1, ncenters(iwann)
            prodw(iat) = 1.0_dp
            do i = 1, ncenters(iwann)
               if(i == iat) cycle
               prodw(iat) = prodw(iat) * distw(i)
            end do
         end do
         do iat = 1, ncenters(iwann)
            charge(Zatoms(iat,iwann)) = charge(Zatoms(iat,iwann)) + 2.0_dp * wannocc(wann_list(iwann)) * prodw(iat)  / (sum(prodw))
         end do
      end do

!      do iat = 1, atoms%nat
!         print *,'Charge of atom ', iat,' is : ', charge(iat)
!      end do
      if(iproc == 0) print *,'Total Wannier charge is : ',sum(charge)
      if(iproc == 0) print *,'Maximal charge is: ',maxval(charge), 'on atom :', maxloc(charge)
      if(iproc == 0) print *,'Minimal charge is: ',minval(charge), 'on atom :', minloc(charge)  

      open(22,file='Wannier_charge.dat',status='unknown')
      do iat = 1, atoms%nat
         write(22, '(E14.6, 2x, E14.6)') charge(iat)/maxval(charge), charge(iat)
!         write(22, '(E14.6, 2x, E14.6)') (charge(iat)-minval(charge))/(maxval(charge)-minval(charge)), charge(iat)
      end do
      close(22)

!Character analysis
      call read_amn_header(seedname,nproj,nband_old,nkpt_old)

      allocate(amn(nband,nproj),stat=i_stat)
      call memocc(i_stat,amn,'amn',subname)
      allocate(tmatrix(nwann,nproj),stat=i_stat)
      call memocc(i_stat,tmatrix,'tmatrix',subname)
      allocate(l(nproj),stat=i_stat)
      call memocc(i_stat,l,'l',subname)
      allocate(mr(nproj),stat=i_stat)
      call memocc(i_stat,mr,'mr',subname)

      call read_amn(seedname,amn,nproj,nband,nkpt_old)
      call read_proj(seedname, nkpt_old, nproj, l, mr)

      call dgemm('N','N',nwann,nproj,nband,1.0d0,umn(1,1),max(1,nwann),&
      &        amn(1,1),max(1,nband),0.0d0,tmatrix(1,1),max(1,nwann))

!DEBUG
!print *,(sum(tmatrix(iat,:)),iat=1,nwann)
!do iat=1,nwann
!   tel = 0.0d0
!   do i_all =1, nproj
!      tel = tel + tmatrix(iat,i_all)**2
!   end do
!   print *,tel
!end do
!DEBUG

      call character_list(nwann,nproj,tmatrix,plotwann,ncenters,wann_list,l,mr) 

      i_all = -product(shape(l))*kind(l)
      deallocate(l,stat=i_stat)
      call memocc(i_stat,i_all,'l',subname)   
      i_all = -product(shape(mr))*kind(mr)
      deallocate(mr,stat=i_stat)
      call memocc(i_stat,i_all,'mr',subname)   
      i_all = -product(shape(amn))*kind(amn)
      deallocate(amn,stat=i_stat)
      call memocc(i_stat,i_all,'amn',subname)   
      i_all = -product(shape(tmatrix))*kind(tmatrix)
      deallocate(tmatrix,stat=i_stat)
      call memocc(i_stat,i_all,'tmatrix',subname)   
      i_all = -product(shape(distw))*kind(distw)
      deallocate(distw,stat=i_stat)
      call memocc(i_stat,i_all,'distw',subname)   
      i_all = -product(shape(charge))*kind(charge)
      deallocate(charge,stat=i_stat)
      call memocc(i_stat,i_all,'charge',subname)   
      i_all = -product(shape(prodw))*kind(prodw)
      deallocate(prodw,stat=i_stat)
      call memocc(i_stat,i_all,'prodw',subname)   

!DEBUG Test distribution on all the atoms
!!     print *,'######################ENTERING DEBUG#################'
!!      allocate(distw(atoms%nat),stat=i_stat)
!!      call memocc(i_stat,distw,'distw',subname)
!!      allocate(charge(atoms%nat),stat=i_stat)
!!      call memocc(i_stat,charge,'charge',subname)
!!      allocate(prodw(atoms%nat),stat=i_stat)
!!      call memocc(i_stat,prodw,'prodw',subname)
!!      allocate(cxyz2(3,nwann),stat=i_stat)
!!      call memocc(i_stat,cxyz2,'cxyz2',subname)
!!      allocate(list(nwann),stat=i_stat)
!!      call memocc(i_stat,list,'list',subname)
!!      do iwann = 1, nwann
!!         list(iwann) = iwann
!!      end do 
!!      call read_centers(iproc,nwann,nwann,atoms%nat,seedname,list,cxyz2,rxyz_wann,.false.,sprd)
!!      charge = 0.0_dp
!!      do iwann = 1, nwann
!!         distw = 0.0_dp
!!         do iat = 1, atoms%nat
!!            distw(iat) = (rxyz_wann(1,iat)-cxyz2(1,iwann))**2 +&
!!                        (rxyz_wann(2,iat)-cxyz2(2,iwann))**2 +&
!!                        (rxyz_wann(3,iat)-cxyz2(3,iwann))**2
!!         end do
!!         prodw = 0.0_dp
!!         do iat = 1, atoms%nat
!!            prodw(iat) = 1.0_dp
!!            do i = 1, atoms%nat
!!               if(i == iat) cycle
!!print *,'distw(i),prodw(iat)',distw(i),prodw(iat)
!!               prodw(iat) = prodw(iat) * distw(i)
!!            end do
!!         end do
!!         do iat = 1, atoms%nat
!!            charge(iat) = charge(iat) + 2.0_dp * wannocc(iwann) * prodw(iat)  / (sum(prodw))
!!         end do
!!      end do
!!      if(iproc == 0) print *,'Total Wannier(all) charge is : ',sum(charge)
!!      if(iproc == 0) print *,'Maximal charge is: ',maxval(charge), 'on atom :', maxloc(charge)
!!      if(iproc == 0) print *,'Minimal charge is: ',minval(charge), 'on atom :', minloc(charge)
!!      open(22,file='Wannier_charge2.dat',status='unknown')
!!      do iat = 1, atoms%nat
!!         write(22, '(E14.6, 2x, E14.6)') charge(iat)/maxval(charge), charge(iat)
!!      end do
!!      close(22)
!!      i_all = -product(shape(distw))*kind(distw)
!!      deallocate(distw,stat=i_stat)
!!      call memocc(i_stat,i_all,'distw',subname)   
!!      i_all = -product(shape(charge))*kind(charge)
!!      deallocate(charge,stat=i_stat)
!!      call memocc(i_stat,i_all,'charge',subname)   
!!      i_all = -product(shape(prodw))*kind(prodw)
!!      deallocate(prodw,stat=i_stat)
!!      call memocc(i_stat,i_all,'prodw',subname)   
!!!END DEBUG

      i_all = -product(shape(wannocc))*kind(wannocc)
      deallocate(wannocc,stat=i_stat)
      call memocc(i_stat,i_all,'wannocc',subname)   

      !DOS of the sprd
      allocate(types(2,nwann))
      types = 1
      call wannier_dos('sprd.dat',1,plotwann,2,types,sprd)
      deallocate(types)

!      call scalar_kmeans_diffIG(0,maxval(sprd(1:plotwann-1))*1.0d-1,plotwann,sprd,'spread',nsprd,buf)
      call scalar_kmeans_diffIG(iproc,0,sprddiff,plotwann,sprd,'spread',nsprd,buf)

      ! Transform rxyz_wann from angstrom to bohr
      do iat = 1, atoms%nat
        do i = 1, 3
           rxyz_wann(i,iat) = rxyz_wann(i,iat) / b2a
        end do
      end do
      ! Transform cxyz from angstrom to bohr
      do iwann = 1, plotwann
        do i = 1, 3
           cxyz(i,iwann) = cxyz(i,iwann) / b2a
        end do
      end do
      ! Transform also the refpos if necessary
      if(units=='angs') then
        do i=1, 3
           refpos(i) = refpos(i) / b2a
        end do
      end if

!DEBUG
!!open(22,file='test.xyz', status='unknown')
!!write(22,*) plotwann+atoms%nat
!!write(22,*)
!!do iat = 1, plotwann
!!   ! Now print the information
!!   write(22,'(A,3(2x,E14.6))'),'X',(cxyz(i,iat), i=1,3)
!!end do
!!do iat = 1, atoms%nat
!!   ! Now print the information
!!   write(22,'(A,3(2x,E14.6))'),'B',(rxyz_wann(i,iat), i=1,3)
!!end do
!!close(22)

! For now choosing the reference point by hand
!!    refpos(1) = rxyz_wann(1,78) ; refpos(2) = rxyz_wann(2,78) ; refpos(3) = rxyz_wann(3,78)   !Pole of 14-6
!    refpos(1) = rxyz_wann(1,64) ; refpos(2) = rxyz_wann(2,64) ; refpos(3) = rxyz_wann(3,64)
!    refpos(1) = rxyz_wann(1,1) ; refpos(2) = rxyz_wann(2,1) ; refpos(3) = rxyz_wann(3,1)
     
   ! Non charged
!   refpos(1) = (rxyz_wann(1,22) + rxyz_wann(1,43) + rxyz_wann(1,39) + rxyz_wann(1,26) + rxyz_wann(1,6)) / 5.0
!   refpos(2) = (rxyz_wann(2,22) + rxyz_wann(2,43) + rxyz_wann(2,39) + rxyz_wann(2,26) + rxyz_wann(2,6)) / 5.0
!   refpos(3) = (rxyz_wann(3,22) + rxyz_wann(3,43) + rxyz_wann(3,39) + rxyz_wann(3,26) + rxyz_wann(3,6)) / 5.0
!   refpos(1) = (rxyz_wann(1,47) + rxyz_wann(1,25) + rxyz_wann(1,8) + rxyz_wann(1,6) + rxyz_wann(1,22) + rxyz_wann(1,50)) / 6.0
!   refpos(2) = (rxyz_wann(2,47) + rxyz_wann(2,25) + rxyz_wann(2,8) + rxyz_wann(2,6) + rxyz_wann(2,22) + rxyz_wann(2,50)) / 6.0
!   refpos(3) = (rxyz_wann(3,47) + rxyz_wann(3,25) + rxyz_wann(3,8) + rxyz_wann(3,6) + rxyz_wann(3,22) + rxyz_wann(3,50)) / 6.0
!END DEBUG

      ! Output the posref for visual check
      open(22,file='pos_ref.xyz', status='unknown')
      write(22,'(I4)') atoms%nat+1
      write(22,*) !skip this line
      write(22,'(A,3(2x,E14.6))') 'X',(refpos(i),i=1,3)
      do i = 1, atoms%nat
         write(22,'(A,3(2x,E14.6))')atoms%atomnames(atoms%iatype(i)),rxyz_wann(1,i),rxyz_wann(2,i),rxyz_wann(3,i)
      end do
      close(22)

      ! Calculate Center of mass
      CM = 0.0_dp
      do iat = 1, atoms%nat
         do j = 1, 3 
            CM(j) = CM(j) + rxyz_wann(j,iat) / real(atoms%nat,kind=8)
         end do
      end do

      !Calculate the radius of the sphere (choose it to be the biggest distance from the CM)
      rad= 0.0_dp
      do iat = 1, atoms%nat
         dist = 0.0_dp
         do j = 1, 3
            dist = dist + (rxyz_wann(j,iat) - CM(j))**2
         end do
         rad = max(dist, rad)
      end do
      rad =sqrt(rad)

      allocate(proj(atoms%nat,3),stat=i_stat)
      call memocc(i_stat,proj,'proj',subname)
      allocate(projC(plotwann,3),stat=i_stat)
      call memocc(i_stat,projC,'projC',subname)
      allocate(nfacets(plotwann))
      allocate(facets(plotwann,maxval(ncenters)*(maxval(ncenters)-1)/2,3))
      allocate(vertex(plotwann,maxval(ncenters)*(maxval(ncenters)-1)/2,3))

      ! Do stereographic projection of atoms and Wannier centers
      call stereographic_projection(atoms%nat,rxyz_wann,refpos, CM, rad, proj, normal, NeglectPoint)
      call stereographic_projection(plotwann,cxyz,refpos, CM, rad, projC, normal, CNeglectPoint)
      call shift_stereographic_projection(plotwann,projC,atoms%nat,proj)
      call write_stereographic_projection(22, 'proj.xyz    ', atoms, proj, NeglectPoint) 

!DEBUG
!!   ! Now open file that will contain the stereographic projection
!!   open(22,file='proj_Wan.xyz', status='unknown')
!!   if (CNeglectPoint .ne. 0 .and. NeglectPoint .ne. 0)then
!!      write(22,*) atoms%nat+plotwann-2
!!   else if((CNeglectPoint .ne. 0 .and. NeglectPoint == 0) .or. (CNeglectPoint == 0 .and. NeglectPoint .ne. 0) ) then
!!      write(22,*) atoms%nat+plotwann-1
!!   else
!!      write(22,*) plotwann+atoms%nat
!!   end if
!!   write(22,*)
!!
!!   do iat = 1, plotwann
!!      ! Now print the information
!!      if(iat == CNeglectPoint) then
!!         write(22,'(A,3(2x,E14.6))'),'#   '//'X',(projC(iat,i), i=1,3)
!!      else
!!         write(22,'(A,3(2x,E14.6))'),'X',(projC(iat,i), i=1,3)
!!      end if
!!   end do
!!   do iat = 1, atoms%nat
!!      ! Now print the information
!!      if(iat == NeglectPoint) then
!!         write(22,'(A,3(2x,E14.6))'),'#   '//'B',(proj(iat,i), i=1,3)
!!      else
!!         write(22,'(A,3(2x,E14.6))'),'B',(proj(iat,i), i=1,3)
!!      end if
!!   end do
!!
!!   close(22)
!END DEBUG

      !Must warn if a Wannier center is on the projection reference
      if(CNeglectPoint .ne. 0) then
         write(*,*) 'The Wannier center ',CNeglectPoint,'is on the refence point of the projection.'
         write(*,*) 'Surfaces will be deformed'
         call mpi_finalize(ierr)
         stop
      end if

      call build_stereographic_graph_facets(atoms%nat,plotwann,4.0d0,rxyz_wann,ncenters,Zatoms,nfacets,facets,vertex)

!DEBUG
!do iwann = 1, plotwann
!  do i = 1, nfacets(iwann)
!    print *, 'iwann, Facets', iwann, facets(iwann,i,:)
!  end do
!end do
!END DEBUG

      call output_stereographic_graph(atoms%nat,proj,projC,plotwann,ncenters,Zatoms,nfacets,facets,vertex,normal,NeglectPoint)

      i_all = -product(shape(sprd))*kind(sprd)
      deallocate(sprd,stat=i_stat)
      call memocc(i_stat,i_all,'sprd',subname)
      i_all = -product(shape(cxyz))*kind(cxyz)
      deallocate(cxyz,stat=i_stat)
      call memocc(i_stat,i_all,'cxyz',subname)
      i_all = -product(shape(Zatoms))*kind(Zatoms)
      deallocate(Zatoms,stat=i_stat)
      call memocc(i_stat,i_all,'Zatoms',subname)
      i_all = -product(shape(rxyz_wann))*kind(rxyz_wann)
      deallocate(rxyz_wann,stat=i_stat)
      call memocc(i_stat,i_all,'rxyz_wann',subname)
      i_all = -product(shape(proj))*kind(proj)
      deallocate(proj,stat=i_stat)
      call memocc(i_stat,i_all,'proj',subname)
      i_all = -product(shape(projC))*kind(projC)
      deallocate(projC,stat=i_stat)
      call memocc(i_stat,i_all,'projC',subname)
      if(.not. hamilAna) then
         i_all = -product(shape(ncenters))*kind(ncenters)
         deallocate(ncenters,stat=i_stat)
         call memocc(i_stat,i_all,'ncenters',subname)
      end if

     if(iproc==0) then
         write(*,*) '!==================================!'
         write(*,*) '!     Bonding Analysis : DONE      !' 
         write(*,*) '!==================================!'
     end if

  end if
  if (hamilAna) then
     if(iproc==0) then
        write(*,*) '!==================================!'
        write(*,*) '!     Hamiltonian Analysis :       !' 
        write(*,*) '!==================================!'
     end if


      call read_nrpts_hamiltonian(iproc,seedname,nrpts) 

      allocate(ham(nrpts,nwann,nwann),stat=i_stat)
      call memocc(i_stat,ham,'ham',subname)

      call read_hamiltonian(iproc,nrpts,nwann,seedname,ham)


      !Eliminate the unoccupied states
      allocate(hamr(nrpts,plotwann,plotwann),stat=i_stat)
      call memocc(i_stat,hamr,'hamr',subname)
      do i = 1, nrpts
        do iiwann = 1, plotwann
              iw1 = wann_list(iiwann)
           do iwann = 1, plotwann
               iw2 = wann_list(iwann)
               hamr(i,iiwann,iwann) = ham(i,iw1,iw2)
           end do
        end do
      end do

     if(.not.bondAna)then
        allocate(buf(plotwann),stat=i_stat)
        call memocc(i_stat,buf,'buf',subname)
        buf = 1
     end if

     ! Diagonal fo the hamiltonian matrix
     allocate(diag(nrpts,nwann),stat=i_stat)
     call memocc(i_stat,diag,'diag',subname)
     do i = 1, nrpts
        do iwann = 1, nwann
           diag(i,iwann) = ham(i,iwann,iwann)
        end do
     end do
     allocate(types(nrpts,nwann),stat=i_stat)
     call memocc(i_stat,types,'types',subname)
     do i = 1, nrpts
        do iwann = 1, nwann
           notocc = .true.
           do iiwann = 1, plotwann
              if(iwann ==  wann_list(iiwann)) then
                 iw1 = buf(iiwann)
                 notocc = .false. 
                 exit
              end if
           end do
           if(notocc) iw1 = nsprd + 1
           types(i,iwann) = iw1
        end do
     end do

     allocate(eigen(1,nband),stat=i_stat)
     call memocc(i_stat,eigen,'eigen',subname)
     call read_eigenvalues(trim(seedname)//'.eig',nband,1,eigen)
     call wannier_projected_dos('Wannier_projected_dos.dat',nrpts,nwann,nband,umn,nsprd+1,types,eigen)
     call wannier_dos('Wannier_dos.dat',nrpts,nwann,nsprd+1,types,diag)
     i_all = -product(shape(eigen))*kind(eigen)
     deallocate(eigen,stat=i_stat)
     call memocc(i_stat,i_all,'eigen',subname)

       if(iproc == 0) write(*,'(A)') 'Diagonal of the Hamiltonian (without empty WFs)'
        allocate(diagT(plotwann),stat=i_stat)
        call memocc(i_stat,diagT,'diagT',subname)
        do i = 1, nrpts
           do iwann = 1, plotwann
              if(bondAna .and. iproc==0) then
                 write(*,'(i4,2x,i4,2x,E14.6,2x,i4)')iwann,iwann,ham(i,iwann,iwann),ncenters(iwann)
              else if(iproc==0) then
                 write(*,'(i4,2x,i4,2x,E14.6)')iwann,iwann,ham(i,iwann,iwann)
              end if
              diagT(iwann) = hamr(i,iwann,iwann)
           end do
        end do

     !Deallocate buf, because it is allocated again in scalar_kmeans_diffIG
     i_all=-product(shape(buf))*kind(buf)
     deallocate(buf,stat=i_stat)
     call memocc(i_stat,i_all,'buf',subname)
     if(.not. bondAna) nsprd = 0  !if we didn't do bonding analysis, find here the best for the hamiltonian
     call scalar_kmeans_diffIG(iproc,nsprd,enediff,plotwann,diagT,'diagonal',ndiag,buf)

!! DEBUG
!! Diagonalizing the Hamiltonian yields the Kohn-Sham eigenvalues and orbitals (U^-1)
!!   allocate(Uham(nwann*(nwann+1)/2),stat=i_stat)
!!   call memocc(i_stat,Uham,'Uham',subname)
!!
!!   ! Copy the upper half of the hamiltonian
!!   !do i = 1, nrpts     ! Right now only one cell
!!   ind = 0
!!   do iiwann = 1, nwann
!!      do iwann = iiwann, nwann
!!           ind = ind +1
!!           Uham(ind) = hamr(1,iiwann,iwann)
!!           if(iproc==0)print *,iiwann,iwann,hamr(1,iiwann,iwann)
!!      end do
!!   end do
!!
!!   allocate(evalues(nwann),stat=i_stat)
!!   call memocc(i_stat,evalues,'evalues',subname)
!!   allocate(evectors(nwann,nwann),stat=i_stat)
!!   call memocc(i_stat,evectors,'evectors',subname)
!!   allocate(work(1+6*nwann+nwann**2),stat=i_stat)
!!   call memocc(i_stat,work,'work',subname)
!!   allocate(iwork(3+5*nwann),stat=i_stat)
!!   call memocc(i_stat,iwork,'iwork',subname)
!!
!!   !Now we can solve the matrix: if(ierr == 0), evectors and evalues contains the eigenvectors and eigenvalues 
!!   call dspevd('V', 'L', nwann, Uham, evalues, evectors, nwann, work, 1+6*nwann+nwann**2, iwork, 3+5*nwann, ierr)
!!
!!   i_all=-product(shape(work))*kind(work)
!!   deallocate(work,stat=i_stat)
!!   call memocc(i_stat,i_all,'work',subname)
!!   i_all=-product(shape(iwork))*kind(iwork)
!!   deallocate(iwork,stat=i_stat)
!!   call memocc(i_stat,i_all,'iwork',subname)
!!
!!
!!   !Print the eigenvalues and eigenvectors
!!   if(iproc == 0) then
!!     do iwann = 1, nwann
!!        write(*,*) 'Eigenvalue:',evalues(iwann)
!!        write(*,*) evectors(iwann,:)
!!     end do
!!   end if
!!
!!   i_all=-product(shape(evalues))*kind(evalues)
!!   deallocate(evalues,stat=i_stat)
!!   call memocc(i_stat,i_all,'evalues',subname)
!!   i_all=-product(shape(evectors))*kind(evectors)
!!   deallocate(evectors,stat=i_stat)
!!   call memocc(i_stat,i_all,'evectors',subname)
!!   i_all=-product(shape(Uham))*kind(Uham)
!!   deallocate(Uham,stat=i_stat)
!!   call memocc(i_stat,i_all,'Uham',subname)
! END DEBUG

      i_all=-product(shape(types))*kind(types)
      deallocate(types,stat=i_stat)
      call memocc(i_stat,i_all,'types',subname)
      i_all=-product(shape(diagT))*kind(diagT)
      deallocate(diagT,stat=i_stat)
      call memocc(i_stat,i_all,'diagT',subname)
      i_all=-product(shape(hamr))*kind(hamr)
      deallocate(hamr,stat=i_stat)
      call memocc(i_stat,i_all,'hamr',subname)
      i_all=-product(shape(ham))*kind(ham)
      deallocate(ham,stat=i_stat)
      call memocc(i_stat,i_all,'ham',subname)
      i_all=-product(shape(diag))*kind(diag)
      deallocate(diag,stat=i_stat)
      call memocc(i_stat,i_all,'diag',subname)
      if(bondAna)then
         i_all = -product(shape(ncenters))*kind(ncenters)
         deallocate(ncenters,stat=i_stat)
         call memocc(i_stat,i_all,'ncenters',subname)
      end if

      if(iproc==0) then 
         write(*,*) '!==================================!'
         write(*,*) '!     Hamiltonian Analysis : DONE  !' 
         write(*,*) '!==================================!'
      end if

  end if
  if(WannCon) then
     !###########################################################################
     ! Set-up number of states used in the Wannier construction
     !###########################################################################
     nvirtu = nband
     nvirtd = 0
     if (input%nspin==2) nvirtd=0!nvirtu
     call orbitals_descriptors(iproc,nproc,nvirtu+nvirtd,nvirtu,nvirtd, &
         & orbs%nspin,orbs%nspinor,orbs%nkpts,orbs%kpts,orbs%kwgts,orbsw)
     call orbitals_communicators(iproc,nproc,Glr,orbsw,commsw)

     nvirtu = n_virt
     nvirtd = 0
     if (input%nspin==2) nvirtd=0!nvirtu
     call orbitals_descriptors(iproc,nproc,nvirtu+nvirtd,nvirtu,nvirtd, &
         & orbs%nspin,orbs%nspinor,orbs%nkpts,orbs%kpts,orbs%kwgts,orbsv)

     !##########################################################################
     ! Read the Wavefunctions
     !##########################################################################

     ! Wavefunctions calculated by BigDFT already are normalized.
     if (iproc==0) then
        write(*,*) '!==================================!'
        write(*,*) '!  Reading the wavefunctions       !'
        write(*,*) '!==================================!'
     end if

     call timing(iproc,'CrtProjectors ','ON')

     ! assign the input_wf_format
     iformat = WF_FORMAT_NONE
     select case (filetype)
     case ("ETSF","etsf")
        iformat = WF_FORMAT_ETSF
     case ("BIN","bin")
        iformat = WF_FORMAT_BINARY
     case default
        if (iproc == 0) write(*,*)' WARNING: Missing specification of wavefunction files'
        stop
     end select

     npsidim=max((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbsw%norbp*orbsw%nspinor,sum(commsw%ncntt(0:nproc-1)))
     allocate(psi(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,orbsw%norbp*orbsw%nspinor),stat=i_stat)
     call memocc(i_stat,psi,'psi',subname)
     allocate(rxyz_old(3,atoms%nat),stat=i_stat)  
     call memocc(i_stat,rxyz_old,'rxyz_old',subname)
     ! For the occupied orbitals, need to modifify norbp,isorb to match the total distributed scheme
     orbs%norbp = n_occ - orbsw%isorb
     if (orbsw%isorb + orbsw%norbp < n_occ ) orbs%norbp = orbsw%norbp
     if(orbsw%isorb > n_occ) orbs%norbp = 0
     orbs%isorb = orbsw%isorb
     if(associated(orbs%iokpt)) then
        i_all = -product(shape(orbs%iokpt))*kind(orbs%iokpt)
        deallocate(orbs%iokpt,stat=i_stat)
        call memocc(i_stat,i_all,'orbs%iokpt',subname)
     end if
     allocate(orbs%iokpt(orbs%norbp),stat=i_stat)
     call memocc(i_stat,orbs%iokpt,'orbs%iokpt',subname)
     orbs%iokpt=1
     if(orbs%norbp > 0) then
        if(associated(orbs%eval)) nullify(orbs%eval)
        allocate(orbs%eval(orbs%norb*orbs%nkpts), stat=i_stat)
        call memocc(i_stat,orbs%eval,'orbs%eval',subname)
        filename=trim(input%dir_output) // 'wavefunction'
        call readmywaves(iproc,filename,iformat,orbs,Glr%d%n1,Glr%d%n2,Glr%d%n3,input%hx,input%hy,input%hz,atoms,rxyz_old,rxyz,  & 
           Glr%wfd,psi(1,1))
        i_all = -product(shape(orbs%eval))*kind(orbs%eval)
        deallocate(orbs%eval,stat=i_stat)
        call memocc(i_stat,i_all,'orbs%eval',subname)
     
     end if

     ! For the non-occupied orbitals, need to change norbp,isorb
     orbsv%norbp = orbsw%isorb + orbsw%norbp - n_occ
     if (orbsw%isorb + orbsw%norbp < n_occ ) orbsv%norbp = 0
     if (orbsw%isorb > n_occ) orbsv%norbp = orbsw%norbp
     orbsv%isorb = 0
     if(orbsw%isorb >= n_occ) orbsv%isorb = orbsw%isorb - n_occ
     if(associated(orbsv%iokpt)) then
        i_all = -product(shape(orbsv%iokpt))*kind(orbsv%iokpt)
        deallocate(orbsv%iokpt,stat=i_stat)
        call memocc(i_stat,i_all,'orbsv%iokpt',subname)
     end if
     allocate(orbsv%iokpt(orbsv%norbp),stat=i_stat)
     call memocc(i_stat,orbsv%iokpt,'orbsv%iokpt',subname)
     orbsv%iokpt=1

     ! read unoccupied wavefunctions
     if(orbsv%norbp > 0) then
        filename=trim(input%dir_output) // 'virtuals'
        if(associated(orbsv%eval)) nullify(orbsv%eval)
        allocate(orbsv%eval(orbsv%norb*orbsv%nkpts), stat=i_stat)
        call memocc(i_stat,orbsv%eval,'orbsv%eval',subname)
        call readmywaves(iproc,filename,iformat,orbsv,Glr%d%n1,Glr%d%n2,Glr%d%n3,input%hx,input%hy,input%hz,atoms,rxyz_old,rxyz,  & 
           Glr%wfd,psi(1,1+orbs%norbp),virt_list)
        i_all = -product(shape(orbsv%eval))*kind(orbsv%eval)
        deallocate(orbsv%eval,stat=i_stat)
        call memocc(i_stat,i_all,'orbsv%eval',subname)
     end if

     i_all = -product(shape(rxyz_old))*kind(rxyz_old)
     deallocate(rxyz_old,stat=i_stat)
     call memocc(i_stat,i_all,'rxyz_old',subname)
     i_all = -product(shape(virt_list))*kind(virt_list)
     deallocate(virt_list,stat=i_stat)
     call memocc(i_stat,i_all,'virt_list',subname)


     call timing(iproc,'CrtProjectors ','OF')

     if (iproc==0) then
        write(*,*) '!===================================!'
        write(*,*) '!  Reading the wavefunctions : DONE !'
        write(*,*) '!===================================!'
     end if

     !#########################################################
     ! Construct the Wannier functions
     !#########################################################

     if(iproc==0) then
        write(*,*) '!==================================!'
        write(*,*) '!     Constructing WFs :           !'
        write(*,*) '!==================================!'
     end if


     allocate(wann(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f),stat=i_stat)
     call memocc(i_stat,wann,'wann',subname) 
     allocate(wannr(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i),stat=i_stat)
     call memocc(i_stat,wannr,'wannr',subname) 
     call initialize_work_arrays_sumrho(Glr,w)


     ! Separate plotwann
!!     allocate(plotwann_par(0:nproc-1),stat=i_stat)
!!     call memocc(i_stat,plotwann_par,'plotwann_par',subname)
!!     allocate(isplotwann_par(0:nproc-1),stat=i_stat)
!!     call memocc(i_stat,isplotwann_par,'isplotwann_par',subname)
!!     call parallel_repartition_with_kpoints(nproc,1,plotwann,plotwann_par)
!!     ntot=0
!!     do jproc=0,nproc-1
!!        isplotwann_par(jproc)=ntot
!!        ntot=ntot+plotwann_par(jproc)
!!     end do
!!     isplotwann = isplotwann_par(iproc) 
!!     plotwannp = plotwann_par(iproc)
!!     i_all = -product(shape(plotwann_par))*kind(plotwann_par)
!!     deallocate(plotwann_par,stat=i_stat)
!!     call memocc(i_stat,i_all,'plotwann_par',subname)
!!     i_all = -product(shape(isplotwann_par))*kind(isplotwann_par)
!!     deallocate(isplotwann_par,stat=i_stat)
!!     call memocc(i_stat,i_all,'isplotwann_par',subname)


     ! Now construct the WFs
     ifile = 12 + iproc
     do iiwann = 1, nwannCon
        iwann = ConstList(iiwann)
        if(iwann == 0) stop 'this should not happen' 
        call to_zero(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f, wann(1))
        do iband = 1, orbsw%norbp
           do i = 1, Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
              wann(i) = wann(i) + umn(iwann,iband+orbsw%isorb) * psi(i,iband)
           end do
        end do
        ! Construction of the Wannier function.
        call mpiallred(wann(1),Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,MPI_SUM,MPI_COMM_WORLD,ierr)
        if(iproc == 0) then
           write(num,'(i4.4)') iwann
           if(outputype == 'cube') then
              !Put it in interpolating scaling functions
              call daub_to_isf(Glr,w,wann(1),wannr)
              call write_wannier_cube(ifile,trim(seedname)//'_'//num//'.cube',atoms,Glr,input,rxyz,wannr)
           else
             if(trim(outputype)=='bin') then
                open(ifile, file=trim(seedname)//'_'//num//'.bin', status='unknown',form='unformatted')
                call writeonewave(ifile,.false.,iiwann,Glr%d%n1,Glr%d%n2,Glr%d%n3,input%hx,input%hy,input%hz,atoms%nat,rxyz,  & 
                   Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keygloc(1,1),Glr%wfd%keyvloc(1),  & 
                   Glr%wfd%nseg_f,Glr%wfd%nvctr_f,Glr%wfd%keygloc(1,Glr%wfd%nseg_c+1),Glr%wfd%keyvloc(Glr%wfd%nseg_c+1), & 
                     wann(1),wann(Glr%wfd%nvctr_c+1), 0.d0)
             else
                stop 'ETSF not implemented yet'                
                ! should be write_wave_etsf  (only one orbital)
                call write_waves_etsf(iproc,trim(seedname)//'_'//num//'.etsf',orbs,Glr%d%n1,Glr%d%n2,Glr%d%n3,&
                     input%hx,input%hy,input%hz,atoms,rxyz,Glr%wfd,wann)
             end if 
           end if
        end if
     end do  ! closing loop on iwann

     if(iproc==0) then
        write(*,*) '!==================================!'
        write(*,*) '!     Constructing WFs : DONE      !'
        write(*,*) '!==================================!'
     end if


     call deallocate_work_arrays_sumrho(w)
     call deallocate_orbs(orbsv,subname)
     call deallocate_orbs(orbsw,subname)
     call deallocate_comms(commsw,subname)


     i_all = -product(shape(wann))*kind(wann)
     deallocate(wann,stat=i_stat)
     call memocc(i_stat,i_all,'wann',subname)
     i_all = -product(shape(wannr))*kind(wannr)
     deallocate(wannr,stat=i_stat)
     call memocc(i_stat,i_all,'wannr',subname)
     i_all = -product(shape(psi))*kind(psi)
     deallocate(psi,stat=i_stat)
     call memocc(i_stat,i_all,'psi',subname)
  end if

  if(.not. WannCon) then
     i_all = -product(shape(virt_list))*kind(virt_list)
     deallocate(virt_list,stat=i_stat)
     call memocc(i_stat,i_all,'virt_list',subname)    
  end if
  i_all = -product(shape(umn))*kind(umn)
  deallocate(umn,stat=i_stat)
  call memocc(i_stat,i_all,'umn',subname)
  if(bondana .or. hamilana) then
     i_all = -product(shape(buf))*kind(buf)
     deallocate(buf,stat=i_stat)
     call memocc(i_stat,i_all,'buf',subname)
  end if
  i_all = -product(shape(wann_list))*kind(wann_list)
  deallocate(wann_list,stat=i_stat)
  call memocc(i_stat,i_all,'wann_list',subname)
  i_all = -product(shape(rxyz))*kind(rxyz)
  deallocate(rxyz,stat=i_stat)
  call memocc(i_stat,i_all,'rxyz',subname)
  call deallocate_lr(Glr,subname)
  call deallocate_orbs(orbs,subname)
  !call deallocate_atoms_scf(atoms,subname)
  call deallocate_atoms(atoms,subname)
  call free_input_variables(input)

  !#########################################################
  ! Ending timing and MPI
  !#########################################################
  call timing(iproc,'             ','RE')

  call cpu_time(tcpu1)
  call system_clock(ncount1,ncount_rate,ncount_max)
  tel=dble(ncount1-ncount0)/dble(ncount_rate)
  if (iproc == 0) &
    write( *,'(1x,a,1x,i4,2(1x,f12.2))') 'CPU time/ELAPSED time for root process ', iproc,tel,tcpu1-tcpu0
   !finalize memory counting
   call memocc(0,0,'count','stop')
 
  ! Barrier suggested by support for titane.ccc.cea.fr, before finalise.
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 
  call MPI_FINALIZE(ierr)

end program Waco

subroutine Waco_input_variables(iproc,filename,nband,nwann,bondAna,Stereo,hamilAna,WannCon,filetype,nwannCon,refpos,units,&
           sprdfact,sprddiff,enediff)
   use module_base
   use module_types
   use module_input
   implicit none
   integer, intent(in) :: iproc
   character(len=*), intent(in) :: filename
   logical, intent(out) :: bondAna, Stereo, hamilAna, WannCon
   character(len=4), intent(out) :: filetype,units
   integer, intent(out) :: nband,nwann,nwannCon 
   real(gp), intent(out) :: sprdfact,sprddiff,enediff
   real(gp), dimension(3), intent(out) :: refpos
   ! Local variable
   logical :: exists

   ! Open the file
   call input_set_file(iproc,.false.,trim(filename),exists,'Waco Parameters')

   ! Read the number of bands and the number of Wannier functions
   call input_var(nband,'1')
   call input_var(nwann,'1', comment='!Number of bands and Wannier functions')

   ! Reading the first three logicals
   call input_var(bondAna,'T')
   call input_var(Stereo,'F')
   call input_var(hamilAna,'T',comment='! Bonding analysis, Stereographic projection, Hamiltonian analysis')

   ! Read the variables for the bonding analysis
   call input_var(sprdfact, '1.64')
   call input_var(sprddiff, '0.1')
   call input_var(enediff, '0.1', comment='! Spread factor for bonding analysis, intervals for k-means of spread and energy ') 

   ! Read the reference position for the stereographic projection
   call input_var(refpos(1),'0.0')
   call input_var(refpos(2),'0.0')
   call input_var(refpos(3),'0.0')
   call input_var(units,'bohr', comment='! Reference position for stereographic projection, units (bohr, angs)')

   ! Check for the wannier construction
   call input_var(WannCon,'F')
   call input_var(filetype,'cube', comment='! Wannier function construction, type of output file (cube, bin or etsf)')
   call input_var(nwannCon,'1',comment='! number of Wannier to construct (if 0 do all) followed by Wannier list on next &
   &     line (optional)')
   
   !Make defaults
   call input_free()

end subroutine Waco_input_variables

subroutine read_input_waco(filename,nwannCon,Constlist)
  use module_base
  implicit none
  character(len=*), intent(in) :: filename
  integer, intent(in) :: nwannCon
  integer, dimension(nwannCon), intent(out) :: Constlist
  ! Local variable
  integer :: i

  open(22, file=trim(filename),status='old')
  ! Skip first lines
  do i=1,6
     read(22,*)
  end do
  ! read the Constlist
  read(22,*) (Constlist(i),i=1,nwannCon)
  close(22)
  
end subroutine read_input_waco

!>
subroutine read_inter_list(iproc,n_virt, virt_list)

   ! This routine reads the list of virtual orbitals needed

   implicit none

   ! I/O variables
   integer, intent(in) :: n_virt,iproc
   integer, dimension(n_virt), intent(out) :: virt_list

   ! Local variables
   integer :: i,j,ierr 

   open(11, file='input.inter', status='old')

   !   write(*,*) '!==================================!'
   !   write(*,*) '!  Reading virtual orbitals list : !'
   !   write(*,*) '!==================================!'

   do i=1,6
      read(11,*,iostat=ierr) ! Skip first lines       
   end do
   read(11,*,iostat=ierr) (virt_list(j), j=1,n_virt)

   if(ierr < 0) then  !reached the end of file and no virt_list, so generate the trivial one
      do j= 1, n_virt
         virt_list(j) = j
      end do
   end if
   close(11)

   if (iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '!Reading virtual orbitals list done!'
      write(*,*) '!==================================!'
      print *
      print *
   end if

END SUBROUTINE read_inter_list


subroutine read_inter_header(iproc,seedname, filetype, n_occ, pre_check, n_virt_tot, n_virt, w_unk, w_sph, w_ang, w_rad)

   ! This routine reads the first lines of a .inter file

   implicit none

   ! I/O variables
   integer, intent(in) :: iproc
   character, intent(out) :: seedname*16, filetype*4
   integer, intent(out) :: n_occ, n_virt, n_virt_tot
   logical, intent(out) :: w_unk, w_sph, w_ang, w_rad, pre_check

   ! Local variables
   character :: char1*1, char2*1, char3*1, char4*1
   logical :: file_exist
   integer :: ierr

   ! Should check if it exists, if not, make a nice output message
   inquire(file="input.inter",exist=file_exist)
   if (.not. file_exist) then
      if(iproc == 0) then
         write(*,'(A)') 'ERROR : Input file, input.inter, not found !'
         write(*,'(A)') 'CORRECTION: Create or give correct input.inter file.'
      end if
      call mpi_finalize(ierr)
      stop
   end if

   OPEN(11, FILE='input.inter', STATUS='OLD')

   if(iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '!   Reading input.inter header :   !'
      write(*,*) '!==================================!'
   end if

   w_unk=.false.
   w_sph=.false.
   w_ang=.false.
   w_rad=.false.

   ! First line
   read(11,*) seedname
   if(iproc==0)write(*,*) 'System studied : ', trim(seedname)

   ! Second line
   read(11,*) filetype
   if(iproc==0)write(*,*) 'file type : ', filetype

   ! Third line
   read(11,*) n_occ
   if(iproc==0)write(*,'(A30,I4)') 'Number of occupied orbitals :', n_occ

   ! Fourth line
   read(11,*) char1, n_virt_tot, n_virt
   if (char1=='T') then
      pre_check=.true.
      if(iproc==0)write(*,*) 'Pre-check before calculating Amnk and Mmnk matrices'
      if(iproc==0)write(*,'(A38,I4)') 'Total number of unnocupied orbitals :', n_virt_tot
   else
      pre_check=.false.
      if(iproc==0)write(*,*) 'Calculation of Amnk and Mmnk matrices'
      if(iproc==0)write(*,'(A39,I4)') 'Number of chosen unnocupied orbitals :', n_virt
   end if

   ! Fifth line
   read(11,*) char1, char2, char3, char4
   if (char1=='T') then
      w_unk=.true.
      if(iproc==0) write(*,*) 'You want to write a UNKp.s file'
   else if (char1 /= 'F') then
      if(iproc==0) write(*,*) 'Wrong value for w_unk'
      STOP
   end if
   if (char2=='T') then
      w_sph=.true.
      if(iproc==0) write(*,*) 'You want to write .cube files for spherical harmonics'
   else if (char2 .ne. 'F') then
      if(iproc==0) write(*,*) 'Wrong value for w_sph'
      STOP
   end if
   if (char3=='T') then
      w_ang=.true.
      if(iproc==0)write(*,*) 'You want to write .cube files for angular parts of the spherical harmonics'
   else if (char3 .ne. 'F') then
      if(iproc==0) write(*,*) 'Wrong value for w_ang'
      STOP
   end if
   if (char4=='T') then
      w_rad=.true.
      if(iproc==0)write(*,*) 'You want to write .cube files for radial parts of the spherical harmonics'
   else if (char4 .ne. 'F') then
      if(iproc==0) write(*,*) 'Wrong value for w_rad'
      STOP
   end if

   CLOSE(11)

   if(iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '! Reading input.inter header done  !'
      write(*,*) '!==================================!'
      print *
      print *
   end if

end subroutine read_inter_header

subroutine read_umn(iproc,nwann,nband,seedname,umn)
   use module_types
   implicit none
   integer, intent(in) :: iproc
   integer, intent(in) :: nwann, nband
   character(len=*),intent(in) :: seedname
   real(gp),dimension(nwann,nband),intent(out) :: umn
   !Local variables
   logical :: file_exist
   integer :: ierr, nwann_umn,nband_umn,iwann,iband,int1,int2

   if(iproc==0) then
      write(*,*) '!==================================!'
      write(*,*) '!     Reading .umn :               !' 
      write(*,*) '!==================================!'
   end if

   inquire(file=trim(seedname)//'.umn',exist=file_exist)
   if (.not. file_exist) then
      if (iproc==0) then
         write(*,'(A,1x,A)') 'ERROR : Input file,',trim(seedname)//'.umn, not found !'
         write(*,'(A)') 'CORRECTION: Create or give correct input file.'
      end if
      call mpi_finalize(ierr)
      stop
   end if

   open(11, file=trim(seedname)//'.umn', status='OLD')
   read(11,*) nwann_umn, nband_umn

   if(nwann_umn .ne. nwann .or. nband_umn .ne. nband) then
     if(iproc == 0) then
       write(*,'(A,I4)') 'ERROR : number of wannier functions in umn,',nwann_umn
       write(*,'(A,I4)') 'not equal number of Wannier functions used:',nwann
       write(*,'(A,I4)') 'ERROR : number of orbitals in umn,',nband_umn
       write(*,'(A,I4)') 'not equal number of orbitals used ',nband
     end if
     call mpi_finalize(ierr)
     stop 
   end if

   do iwann = 1,nwann
      do iband = 1,nband
         read(11,*) int1, int2, umn(iwann,iband)
      end do
   end do
   close(11)

   if(iproc==0) then
      write(*,*) 'Number of Wannier functions : ',nwann_umn
      write(*,*) 'Number of orbitals used     : ',nband_umn
      write(*,*) '!==================================!'
      write(*,*) '!     Reading .umn : DONE          !' 
      write(*,*) '!==================================!'
   end if
end subroutine read_umn

subroutine read_centers(iproc,nwann,plotwann,natom,seedname,wann_list,cxyz,rxyz_wann,readsprd,sprd)
   use module_types
   implicit none
   logical, intent(in) :: readsprd
   integer, intent(in) :: iproc,plotwann,nwann,natom
   character(len=*),intent(in) :: seedname
   integer, dimension(plotwann),intent(in) :: wann_list
   real(gp),dimension(3,plotwann),intent(out) :: cxyz
   real(gp),dimension(3,natom),intent(out) :: rxyz_wann
   real(gp),dimension(plotwann+1), intent(out) :: sprd
   !Local variables
   character(len=4) :: char1
   logical :: file_exist, commented
   integer :: i,iwann,iiwann

   !read the _centres.xyz file
   inquire(file=trim(seedname)//'_centres.xyz',exist=file_exist)
   if (.not. file_exist) then
      if (iproc==0) then
         write(*,'(A,1x,A)') 'ERROR : Input file,',trim(seedname)//'_centres.xyz, not found !'
         write(*,'(A)') 'CORRECTION: Create or give correct input file.'
      end if
      call mpi_finalize(i)
      stop
   end if

   open(11, file=trim(seedname)//'_centres.xyz', status='OLD')

   !skip first two lines
   read(11,*)
   read(11,*)

   !now read the centers
   iiwann = 0
   do iwann = 1, nwann
      commented = .true.
      do i = 1, plotwann
         if(iwann == wann_list(i)) commented = .false.
      end do
      if (.not. commented) then
         iiwann = iiwann + 1
         read(11,*) char1, cxyz(1,iiwann), cxyz(2,iiwann), cxyz(3,iiwann)
      end if
      if (commented) read(11,*)  !just skip this line
   end do

   !now read atomic positions (already in rxyz but beware shifts...)
   do i = 1, natom
      read(11,*) char1, rxyz_wann(1,i), rxyz_wann(2,i), rxyz_wann(3,i)
   end do

   ! now read the spreads (the additionnal last value is total sprd)
   if(readsprd) then
   print *,'plotwann',plotwann
      do iiwann = 1, plotwann+1
         read(11, *) char1, sprd(iiwann)
      end do
   end if
   close(11)
end subroutine read_centers

subroutine read_nrpts_hamiltonian(iproc,seedname,nrpts)
   implicit none
   integer, intent(in) :: iproc
   character(len=*), intent(in) :: seedname
   integer, intent(out) :: nrpts
   ! Local variables
   integer :: ierr
   logical :: file_exist

   !read the _hr.dat file
   inquire(file=trim(seedname)//'_hr.dat',exist=file_exist)
   if (.not. file_exist) then
      if (iproc==0) then
         write(*,'(A,1x,A)') 'ERROR : Input file,',trim(seedname)//'_hr.dat, not found !'
         write(*,'(A)') 'CORRECTION: Create or give correct input file.'
      end if
      call mpi_finalize(ierr)
      stop
   end if
   
   !read the hamiltonian
   open(11, file=trim(seedname)//'_hr.dat', status='OLD')
   
   !skip first lines
   read(11,*)           !comment line containing the date 
   read(11,*)           !this line contains the number of wannier functions (already known)
   read(11,*) nrpts     ! this line contains the number of Wigner-Seitz grid points
   close(11)

end subroutine read_nrpts_hamiltonian

subroutine read_hamiltonian(iproc,nrpts,nwann,seedname,ham)
   use module_types
   implicit none
   integer, intent(in) :: iproc
   integer, intent(in) :: nrpts
   integer, intent(in) :: nwann
   character(len=*), intent(in) :: seedname
   real(gp), dimension(nrpts,nwann,nwann), intent(out) :: ham
   ! Local variables
   integer :: i,ierr,iwann, iiwann, int1, int2, int3, iw1, iw2
   logical :: file_exist

   !read the _hr.dat file
   inquire(file=trim(seedname)//'_hr.dat',exist=file_exist)
   if (.not. file_exist) then
      if (iproc==0) then
         write(*,'(A,1x,A)') 'ERROR : Input file,',trim(seedname)//'_hr.dat, not found !'
         write(*,'(A)') 'CORRECTION: Create or give correct input file.'
      end if
      call mpi_finalize(ierr)
      stop
   end if
   
   !read the hamiltonian
   open(11, file=trim(seedname)//'_hr.dat', status='OLD')
   
   !skip first lines
   read(11,*)           !comment line containing the date 
   read(11,*)           !this line contains the number of wannier functions (already known)
   read(11,*)           ! this line contains the number of Wigner-Seitz grid points
   do i=1,nrpts,15
      read(11,*)        !skip the block of integers giving the degeneracy of the Wigner-Seitz grid points (15 integers per line)
   end do

   ! Read hamiltonian
   do i = 1, nrpts 
      do iiwann = 1, nwann
         do iwann = 1, nwann
            read(11,*) int1, int2, int3, iw1, iw2, ham(i,iw1,iw2)
         end do
      end do
   end do
   close(11)
   
end subroutine read_hamiltonian

subroutine write_wannier_cube(ifile,filename,atoms,Glr,input,rxyz,wannr)
   use module_types
   implicit none
   character(len=*), intent(in) :: filename
   integer, intent(in) :: ifile
   type(atoms_data),intent(in) :: atoms
   type(locreg_descriptors), intent(in) :: Glr
   type(input_variables),intent(in) :: input
   real(gp),dimension(3,atoms%nat),intent(in) :: rxyz
   real(gp),dimension(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i),intent(in) :: wannr
   ! Local variables
   logical :: perx, pery, perz
   integer :: nbl1,nbr1,nbl2,nbr2,nbl3,nbr3,rem
   integer :: i,j,ix,iy,iz,ind
   
   perx=(Glr%geocode /= 'F')
   pery=(Glr%geocode == 'P')
   perz=(Glr%geocode /= 'F')
   call ext_buffers(perx,nbl1,nbr1)
   call ext_buffers(pery,nbl2,nbr2)
   call ext_buffers(perz,nbl3,nbr3)
   if(nbr1 > 0) nbr1 = nbr1 + 2
   if(nbr2 > 0) nbr2 = nbr2 + 2
   if(nbr3 > 0) nbr3 = nbr3 + 2
   ! Volumetric data in batches of 6 values per line, 'z'-direction first.
   rem=Glr%d%n3i-floor(real(Glr%d%n3i/6))*6
   
   open(ifile, file=filename, status='unknown')
   write(ifile,*) ' CUBE file for ISF field'
   write(ifile,*) ' Case for'
   write(ifile,'(I4,1X,F12.6,2(1X,F12.6))') atoms%nat, real(0.d0), real(0.d0), real(0.d0)
   write(ifile,'(I4,1X,F12.6,2(1X,F12.6))') Glr%d%n1i-(nbl1+nbr1), 0.5_dp*input%hx, real(0.d0),  real(0.d0)
   write(ifile,'(I4,1X,F12.6,2(1X,F12.6))') Glr%d%n2i-(nbl2+nbr2), real(0.d0),  0.5_dp*input%hy, real(0.d0)
   write(ifile,'(I4,1X,F12.6,2(1X,F12.6))') Glr%d%n3i-(nbl3+nbr3), real(0.d0),  real(0.d0),  0.5_dp*input%hz
   do i=1, atoms%nat
      write(ifile,'(I4,1X,F12.6,3(1X,F12.6))') atoms%nzatom(atoms%iatype(i)), real(0.d0), (real(rxyz(j,i)), j=1,3)
   end do
   !do ix=Glr%d%n1i-nbr1,1+nbl1,-1
   !   do iy=Glr%d%n2i-nbr2,1+nbl2,-1
   !      do iz=Glr%d%n3i-nbr3,1+nbl3,-1
   do ix=1+nbl1,Glr%d%n1i-nbr1
      do iy=1+nbl2,Glr%d%n2i-nbr2
         do iz=1+nbl3,Glr%d%n3i-nbr3
            ind = iz*(Glr%d%n2i*Glr%d%n1i) + iy*Glr%d%n1i + ix
            write(ifile,'(E14.6)',advance='no') real(wannr(ind))
            if ( ( (mod(iz+5-rem,6) .eq. 0) .and. (iz .ne. Glr%d%n3i) ) .or. (iz .eq. 1) ) then
               write(ifile,'(a)') ''
            end if
         end do
      end do
   end do
   close(ifile)

end subroutine write_wannier_cube

subroutine scalar_kmeans_diffIG(iproc,nIG,crit,nel,vect,string,nbuf,buf)
  use BigDFT_API
  use module_interfaces
  implicit none
  integer, intent(in) :: nel,nIG,iproc
  real(kind=8),intent(in) :: crit
  real(kind=8), dimension(nel),intent(in) :: vect
  character(len=*),intent(in) :: string
  integer, intent(out) :: nbuf
  integer, dimension(:), pointer :: buf
  ! local variables
  character(len=*), parameter :: subname='scalar_kmeans_diffIG'
  integer :: i_stat, i_all, i, iel, iiel, iter
  real :: r
  real(kind=8) :: minold, maxold
  logical :: same
  integer, allocatable :: degen(:), oldbuf(:)
  real(kind=8), allocatable :: buffers(:), means(:), diff(:)
  

  ! Implement k-means for the diagonal 
  !1) Initial guess
  if(nIG == 0) then ! Maximum Difference between states (determines the number of clusters)
     allocate(buffers(nel),stat=i_stat)
     call memocc(i_stat,buffers,'buffers',subname)
     nbuf = 0
     buffers = 9999.99999_dp
     do iel = 1, nel
        same = .false.
        do iiel = 1, nel
           if(abs(vect(iel) - buffers(iiel)) < crit) then
             same = .true.
           end if
        end do
        if(.not. same) then
           nbuf = nbuf +1
           buffers(nbuf) = vect(iel)
        end if
     end do
  else if (nIG > 0) then  ! Random IG
     nbuf = nIG
     allocate(buffers(nbuf),stat=i_stat)
     call memocc(i_stat,buffers,'buffers',subname)
     do i= 1, nbuf
        call init_random_seed(i)
        call random_number(r)          !returns a random number between [0,1]
        iel = int(r*real(nel))
        if(iel == 0) iel=1
        buffers(i) = vect(iel)
     end do
  end if

  ! Store the initial guess in means
  allocate(means(nbuf),stat=i_stat)
  call memocc(i_stat,means,'means',subname)
  do i = 1, nbuf
     means(i) = buffers(i)
  end do

  allocate(diff(nbuf),stat=i_stat)
  call memocc(i_stat,diff,'diff',subname)
  allocate(degen(nbuf),stat=i_stat)
  call memocc(i_stat,degen,'degen',subname)
  allocate(buf(nel),stat=i_stat)
  call memocc(i_stat,buf,'buf',subname)
  allocate(oldbuf(nel),stat=i_stat)
  call memocc(i_stat,oldbuf,'oldbuf',subname)
 
  buf = 0 
  iter = 0
  loop_iter: do
     iter = iter + 1
     ! Copy the buffers
     do iel = 1, nel
        oldbuf(iel) = buf(iel)
     end do

     !2) Associate the numbers
     do iel = 1, nel
        do i = 1, nbuf
           diff(i) = abs(vect(iel) - means(i))
        end do
        buf(iel) =  minloc(diff,dim=1)
     end do
   
     !4) Stop iterating if assignement does not change (ordering can change?)
     same = .true.
     do iel = 1, nel
        if (buf(iel) .ne. oldbuf(iel)) then
           same = .false.
           exit
        end if
     end do
     if(same) exit loop_iter
 
     !3) Calculate new means
     do i = 1, nbuf
        means(i) = 0.0_dp
        degen(i) = 0
        do iel = 1, nel
           if(buf(iel) .ne. i) cycle
           means(i) = means(i) + vect(iel)
           degen(i) = degen(i) + 1
        end do
          means(i) = means(i)/real(max(degen(i),1),kind=8)
     end do
  end do loop_iter

  if(iproc == 0) then
     write(*,'(A,1x,i4,1x,A)') 'Convergence reached in',iter,'iterations.'
     write(*,'(A,A,A,1x,i4,1x,A)') 'The ',trim(string),' can be clustered in',nbuf,'elements:'
     do i = 1, nbuf
        minold = huge(minold)
        maxold = -huge(maxold)
        do iel = 1, nel
           if(buf(iel) .ne. i) cycle
           if(vect(iel) < minold) minold = vect(iel)
           if(vect(iel) > maxold) maxold = vect(iel)
        end do
        write(*,'(E14.6,3x,A,i4,2x,A,2x,E14.6,2x,A,2x,E14.6)')   means(i), 'with',degen(i),'elements from:',minold,'to',maxold
     end do
  end if

  i_all = -product(shape(means))*kind(means)
  deallocate(means,stat=i_stat)
  call memocc(i_stat,i_all,'means',subname)
  i_all = -product(shape(buffers))*kind(buffers) 
  deallocate(buffers,stat=i_stat)
  call memocc(i_stat,i_all,'buffers',subname)
  i_all = -product(shape(diff))*kind(diff) 
  deallocate(diff,stat=i_stat)
  call memocc(i_stat,i_all,'diff',subname)
  i_all = -product(shape(degen))*kind(degen)
  deallocate(degen,stat=i_stat)
  call memocc(i_stat,i_all,'degen',subname)
  i_all = -product(shape(oldbuf))*kind(oldbuf)
  deallocate(oldbuf,stat=i_stat)
  call memocc(i_stat,i_all,'oldbuf',subname)

end subroutine scalar_kmeans_diffIG

subroutine init_random_seed(shuffler)
  implicit none
  integer, intent(in) :: shuffler
  integer :: i, n, clock
  integer, dimension(:), allocatable :: seed
          
  call random_seed(size = n)
  allocate(seed(n))
  
  call system_clock(count=clock)
  
  seed = clock*shuffler + 37 * (/ (i - 1, i = 1, n) /)
  call random_seed(put = seed)
  
  deallocate(seed)
end subroutine init_random_seed

subroutine stereographic_projection(natom, rxyz, refpos, CM, rad, proj, normal, dcp)
   use BigDFT_API
   use Poisson_Solver
   use module_interfaces
   implicit none
   integer, intent(in) :: natom
   real(gp), dimension(3,natom), intent(in) :: rxyz
   real(gp), dimension(3), intent(inout) :: refpos, CM
   real(gp), intent(in) :: rad
   integer, intent(out) :: dcp                            ! neglected point (atom corresponding to reference)
   real(gp), dimension(natom,3), intent(out) :: proj ! atom positions in the projection
   real(gp), dimension(3), intent(out) :: normal   ! normal to the plane of projection
   !Local variables
   integer :: iat,j,i_stat,i_all
   real(kind=8) :: norm, norm2, theta, dotprod, distsph
   real(kind=8), dimension(:,:), allocatable :: pos
   real(kind=8), dimension(3) :: r,W,NQ,Q,vect,check
   real(kind=8) :: fact,dist
   character(len=*), parameter :: subname='stereographic_projection'
 
 !  write (*,'(A,3(2x,E14.6))')'Center of mass:',(CM(j),j=1,3)
 !  write (*,'(A,2x,E14.6)')'Radius of containing sphere:', rad 

   allocate(pos(natom,3),stat=i_stat)
   call memocc(i_stat,pos,'pos',subname)
 
   ! Now rad is the radius of the sphere
   ! Must displace all atoms to be on the corresponding sphere
   do iat = 1 , natom
      dist = 0.0_dp
      do j = 1, 3
         Q(j) = rxyz(j,iat) - CM(j)
         dist = dist + (rxyz(j,iat) - CM(j))**2
      end do
      dist  = sqrt(dist)
      do j=1, 3
         NQ(j) = Q(j) / dist
         vect(j) =  NQ(j) * rad - Q(j)
         pos(iat,j) = rxyz(j,iat) + vect(j)
      end do
   end do
   
   ! put the reference on the sphere
   dist =  0.0_dp 
   do j = 1, 3
      Q(j) = refpos(j) - CM(j)
      dist = dist + (refpos(j) - CM(j))**2
   end do
   dist  = sqrt(dist)
   do j=1, 3
      NQ(j) = Q(j) / dist
      vect(j) =  NQ(j) * rad - Q(j)
      refpos(j) = refpos(j) + vect(j)
   end do

!DEBUG
!   open(22,file='pos_sph.xyz', status='unknown')
!   write(22,'(I4)') natom+1
!   write(22,*) !skip this line
!   write(22,'(A,3(2x,E14.6))'), 'X',(refpos(j),j=1,3)
!   do iat = 1, natom
!      write(22,'(A,3(2x,E14.6))'),'B',pos(iat,1),pos(iat,2),pos(iat,3) 
!   end do
!   close(22)
!END DEBUG

   !Calculate the vector perpendicular to the plane of projection
   r(1) = refpos(1) - CM(1)
   r(2) = refpos(2) - CM(2)
   r(3) = refpos(3) - CM(3)
  
!   write (*,'(A,3(2x,E14.6))') 'Normal to the plane:',(r(j)/sqrt(r(1)**2+r(2)**2+r(3)**2),j=1,3)
!   write (*,'(A,3(2x,E14.6))') 'Reference position:',refpos

   ! check that reference point is not on an atom
   ! If so, neglect that atom
   dcp = 0
   do iat = 1, natom
      dotprod = 0.0_dp
      norm  = 0.0_dp
      norm2 = 0.0_dp
      do j= 1, 3
         !calculate vector from ref to the points
         W(j) = refpos(j) - pos(iat,j)
         !calculate vector from ref to the CM
         NQ(j) = refpos(j) - CM(j)
         !calculate the vector from CM to points
         Q(j) = pos(iat,j) - CM(j)
         !calculate dot product
         dotprod = dotprod + NQ(j)*W(j)
         norm  = norm  + NQ(j)*NQ(j)
         norm2 = norm2 + W(j)*W(j)
      end do
      norm = sqrt(norm)
      norm2 = sqrt(norm2)
      if(norm2 < 1.0d-10) norm2 = 1.0d-10

      if(norm2 < 1.0d-1 .and. abs(dotprod/(norm*norm2)) < 1.0d-2) then

        if(dcp .ne. 0) stop 'This should not happen'
        dcp = iat
      end if
   end do

   ! Major loop structure that performs the projection
   check = 0.0_dp
   do iat = 1, natom
      dotprod = 0.0_dp
      norm  = 0.0_dp
      norm2 = 0.0_dp
      do j= 1, 3
         !calculate vector from ref to the points
         W(j) = refpos(j) - pos(iat,j)
         !calculate vector from ref to the points
         NQ(j) = refpos(j) - CM(j)
         !calculate the vector from CM to points
         Q(j) = pos(iat,j) - CM(j)
         !calculate dot product
         dotprod = dotprod + NQ(j)*W(j)
         norm  = norm  + NQ(j)*NQ(j)
         norm2 = norm2 + W(j)*W(j)
      end do
      norm = sqrt(norm)
      norm2 = sqrt(norm2)
 
      !Calculate angle
      if(iat == dcp) then
         theta = 0.0_dp
      else
         theta = acos(dotprod / (norm*norm2))
      end if
      distsph = rad * (theta + (theta**3)/3.0_dp)       !this is the new distance used in the contraction
                                               !The real stereographic projection uses tan(theta) = theta + 1/3*theta**3 + ...
   
      !Calculate factor
      if(dot_product(W,W) > 1.0d-10)then
         fact = (Q(1)*r(1) + Q(2)*r(2) + Q(3)*r(3)) / (W(1)*r(1) + W(2)*r(2) + W(3)*r(3))
      else
         fact = 1.0d+10  !the reference and the atom are the same, so this should be at infinity
      end if

      !calculate point on the plane
      do j= 1, 3
         proj(iat,j) = pos(iat,j) - fact * W(j)
      end do
   
      !Now compactify the Riemann sphere
      !Use the distance on the sphere has the distance from the center of graph
      dist=0.0_dp
      do j=1, 3
         dist = dist + (proj(iat,j) - CM(j))**2
      end do
      dist = sqrt(dist)
      if(dist < 1.0d-6) dist = 1.0d-6
      do j=1, 3
         proj(iat,j) = (proj(iat,j)-CM(j)) / dist
         proj(iat,j) = proj(iat,j) * distsph
         check(j) = check(j) + proj(iat,j)
      end do
   end do

!!   do iat = 1, natom
!!      do j= 1, 3
!!         ! Now shift the positions to be contained in the first quadrant
!!         proj2(iat,j) = proj(iat,j)! - minval(proj(:,j))
!!      end do    
!!   end do


   !Normalize the normal to the plane
   do j = 1, 3
      normal(j) = r(j)/sqrt(r(1)**2+r(2)**2+r(3)**2)
   end do

   i_all = -product(shape(pos))*kind(pos)
   deallocate(pos,stat=i_stat)
   call memocc(i_stat,i_all,'pos',subname)

end subroutine stereographic_projection

subroutine write_stereographic_projection(unitnumb, filename, atoms, proj, NeglectPoint)
   use module_types
   implicit none
   integer, intent(in) :: unitnumb
   character(len=12), intent(in) :: filename
   type(atoms_data), intent(in) :: atoms
   integer, intent(in) :: NeglectPoint                   ! neglected point (atom corresponding to reference)
   real(gp), dimension(atoms%nat,3), intent(in) :: proj  ! atom positions in the projection
   !Local variable
   integer :: i, iat

   ! Now open file that will contain the stereographic projection
   open(unitnumb,file=trim(filename), status='unknown')
   if(NeglectPoint .ne. 0) then
      write(unitnumb,*) atoms%nat-1
   else
      write(unitnumb,*) atoms%nat
   end if
   write(unitnumb,*)

   do iat = 1, atoms%nat
      ! Now print the information
      if(iat == NeglectPoint) then
         write(unitnumb,'(A,3(2x,E14.6))')'#   '//atoms%atomnames(atoms%iatype(iat)),(proj(iat,i), i=1,3)
      else
         write(unitnumb,'(A,3(2x,E14.6))')atoms%atomnames(atoms%iatype(iat)),(proj(iat,i), i=1,3)
      end if
   end do
   close(unitnumb)

end subroutine write_stereographic_projection

subroutine shift_stereographic_projection(nwann,wannr,natom,rxyz)
   use module_types
   implicit none

integer, intent(in) :: nwann, natom
real(gp), dimension(natom,3), intent(inout) :: rxyz
real(gp), dimension(nwann,3), intent(inout) :: wannr
! Local variables
integer :: i,j
real(gp), dimension(3) :: shift

shift(1) = minval(rxyz(:,1))
shift(2) = minval(rxyz(:,2))
shift(3) = minval(rxyz(:,3))

! Shift both contributions to be inside the first quadrant
do i = 1, natom
   do j= 1, 3
      rxyz(i,j) = rxyz(i,j) - shift(j)
   end do
end do

do i = 1, nwann
   do j= 1, 3
      wannr(i,j) = wannr(i,j) - shift(j)
   end do
end do

end subroutine shift_stereographic_projection

subroutine build_stereographic_graph_facets(natoms,nsurf,maxbond,rxyz,ncenters,Zatoms,nfacets,facets,vertex)
   use module_interfaces
   use module_types
   implicit none
   integer, intent(in) :: natoms, nsurf
   real(gp), intent(in) :: maxbond
   real(gp), dimension(3,natoms), intent(in) :: rxyz
   integer, dimension(nsurf), intent(in) :: ncenters
   integer, dimension(natoms,nsurf),intent(in) :: Zatoms
   integer, dimension(nsurf), intent(out) :: nfacets
   integer, dimension(nsurf,maxval(ncenters)*(maxval(ncenters) - 1)/2, 3),intent(out) :: facets
   integer, dimension(nsurf,maxval(ncenters)*(maxval(ncenters) - 1)/2, 3),intent(out) :: vertex    !exactly like facets, but in reference to only the surface
   ! Local variables
   integer :: i, j, k, isurf
   real(gp) :: dist

   ! To do so, choose all the triangles formed from the Wannier center and 2 atoms sharing a bond (inside a maximal bond length)
   facets = 0
   do isurf = 1, nsurf

      ! Build the facets
      nfacets(isurf) = 0
      do i = 1, ncenters(isurf)
         do j = i+1, ncenters(isurf)
            dist = 0.0_dp
            do k = 1, 3
               dist = dist + (rxyz(k,Zatoms(i,isurf)) - rxyz(k,Zatoms(j,isurf)))**2
            end do
            if(sqrt(dist) > maxbond) cycle
            nfacets(isurf) = nfacets(isurf) + 1
            facets(isurf,nfacets(isurf),1) = natoms + 1  !corresponds to the Wannier center
            facets(isurf,nfacets(isurf),2) = Zatoms(i,isurf)
            facets(isurf,nfacets(isurf),3) = Zatoms(j,isurf)
            vertex(isurf,nfacets(isurf),1) = ncenters(isurf) + 1  !corresponds to the Wannier center
            vertex(isurf,nfacets(isurf),2) = i
            vertex(isurf,nfacets(isurf),3) = j
         end do
      end do
   end do

end subroutine build_stereographic_graph_facets

subroutine output_stereographic_graph(natoms,proj,projC,nsurf,ncenters,Zatoms,npoly,poly,vertex,normal,NeglectPoint) 
   use module_types
   implicit none
   integer, intent(in) :: natoms, nsurf, NeglectPoint
   real(gp), dimension(natoms,3), intent(in) :: proj                   ! atom position in the projection
   real(gp), dimension(nsurf, 3), intent(in) :: projC                  ! Wannier centers in the projection
   integer, dimension(nsurf), intent(in) :: ncenters, npoly
   integer, dimension(natoms,nsurf),intent(in) :: Zatoms               ! indexes of all the atoms spanned by the Wannier function
   integer, dimension(nsurf,maxval(ncenters)*(maxval(ncenters)-1)/2,3),intent(in) :: poly
   integer, dimension(nsurf,maxval(ncenters)*(maxval(ncenters)-1)/2,3),intent(in) :: vertex
   real(gp), dimension(3), intent(in) :: normal
   ! Local variables
   integer :: isurf,ipts,i, nsurftot,npts,ndecimal,ncent, num_poly, num_poly_tot
   character(len=14) :: surfname
   character(len=20) :: forma
   logical :: condition

   open(23,file='proj.surf', status='unknown')

   !First line is just a comment 
   write(23,'(A)') 'Stereographic graph'

   !Write the dimensions of the box (simple cubic)
   write(23,'(E14.6, 2x, E14.6, 2x, E14.6)') maxval(proj(:,1))-minval(proj(:,1)), 0.0, maxval(proj(:,2))-minval(proj(:,2))  !dxx dyx dyy
   write(23,'(E14.6, 2x, E14.6, 2x, E14.6)') 0.0, 0.0,  maxval(proj(:,3))-minval(proj(:,3))                                   !dzx dzy dzz

   nsurftot = 0
   npts = 0
   num_poly_tot = 0
   do isurf = 1, nsurf
      num_poly = 0
      !MUST eliminate the wanniers that are less then 3 centers
      if(ncenters(isurf) < 3) cycle
      ! or the 3 centers containing the neglected point
      condition = (Zatoms(1,isurf) == NeglectPoint .or. Zatoms(2,isurf) == NeglectPoint .or. Zatoms(3,isurf) == NeglectPoint )
      if(ncenters(isurf) == 3 .and. condition) cycle
      ! Must eliminate the neglected point from other surfaces
      do ipts = 1, npoly(isurf)
         condition = .false.
         do i = 1, 3
            condition = condition .or. poly(isurf,ipts,i) == NeglectPoint
         end do
         if(.not. condition) num_poly = num_poly + 1
         if(.not. condition) num_poly_tot = num_poly_tot + 1
      end do
      if(num_poly > 0)nsurftot = nsurftot + 1
      if(num_poly > 0)npts = npts + ncenters(isurf) + 1
   end do

   !Fourth line: number of surfaces, total num_polys, total num_points
   write(23,'(I4, 2x, I4, 2x, I4)') nsurftot, num_poly_tot, npts

   do isurf = 1, nsurf
      !MUST eliminate the wanniers that are less then 3 centers
      if(ncenters(isurf) < 3) cycle
      ! or the 3 centers containing the neglected point
      condition = (Zatoms(1,isurf) == NeglectPoint .or. Zatoms(2,isurf) == NeglectPoint .or. Zatoms(3,isurf) == NeglectPoint)
      if(ncenters(isurf) == 3 .and. condition) cycle

      ! Determining the format (number of integers) for surface name 
      ndecimal = 1
      ncent = ncenters(isurf)
      do
        if(int(ncent / 10) == 0) exit
        ncent = ncent / 10
        ndecimal = ndecimal + 1
      end do
      write(forma,'(I1)') ndecimal
      forma = '(I'//trim(forma)//')'
      !Name of the surface
      write(surfname,forma) ncenters(isurf)
      surfname = '2e - '//trim(surfname)//'centers'

      !num_polys and num_points (Must eliminate the neglected point from other surfaces)
      condition = .false.
      num_poly = 0
      do ipts = 1, npoly(isurf)
         condition = .false.
         do i = 1, 3
            condition = condition .or. poly(isurf,ipts,i) == NeglectPoint
         end do
         if(.not. condition)num_poly = num_poly + 1
      end do

      if(num_poly > 0) then
         !Write the surface only if there is some polygones
         write(23,'(A)') trim(surfname)

         if(condition) then
            write(23,'(I4, 2x, I4)') num_poly, ncenters(isurf) + 1
            !number of vertices, i_1 i_2 i_3 ... i_n (index of the vertices)
            do ipts = 1, npoly(isurf)
               condition = .false.
               do i = 1, 3
                  condition = condition .or. poly(isurf,ipts,i) == NeglectPoint
               end do
               if(condition) cycle
               write(forma,'(I4)') ncenters(isurf)
               forma = '(I4,'//trim(forma)//'(2x,I4))'
               write(23,trim(forma)) 3, (vertex(isurf,ipts,i),i=1,3) ! Only triangles
            end do
         else
            write(23,'(I4, 2x, I4)') num_poly, ncenters(isurf) + 1
            !number of vertices, i_1 i_2 i_3 ... i_n (index of the vertices)
            do ipts = 1, npoly(isurf)
               condition = .false.
               do i = 1, 3
                  condition = condition .or. poly(isurf,ipts,i) == NeglectPoint
               end do
               if(condition) cycle
               write(forma,'(I4)') ncenters(isurf) + 1
               forma = '(I4,'//trim(forma)//'(2x,I4))'
               write(23,trim(forma)) 3, (vertex(isurf,ipts,i),i=1,3) ! Only triangles
            end do
         end if

         do ipts = 1, ncenters(isurf)
         !   if(Zatoms(ipts,isurf) == NeglectPoint) cycle
            !coordinates of the vertices (x y z) and the normal to the surface at these vertices (nx ny nz)
            write(23,'(5(E14.6, 2x),E14.6)') proj(Zatoms(ipts,isurf),1), proj(Zatoms(ipts,isurf),2), proj(Zatoms(ipts,isurf),3),&
                 (normal(i), i=1,3)
         end do
         write(23,'(5(E14.6, 2x),E14.6)') projC(isurf,1), projC(isurf,2), projC(isurf,3),&
                 (normal(i), i=1,3)
      end if !on num_poly
   end do
   close(23)

end subroutine output_stereographic_graph

subroutine wannier_dos(filename,nrpts,nwann,ntypes,types,diag)
use module_types
implicit none
character(len=*), intent(in) :: filename
integer, intent(in) :: nrpts                                               ! Number of unit cells (r-points)
integer, intent(in) :: nwann                                               ! Number of Wannier functions 
integer, intent(in) :: ntypes                                              ! Number of types of Wannier functions
integer, dimension(nrpts,nwann), intent(in) :: types                       ! Types of the Wannier functions
real(gp), dimension(nrpts,nwann),intent(in) :: diag                        ! Diagonal elements of the Hamiltonian in Wannier basis
!Local variables
integer, parameter :: ndos = 1000
real(gp), parameter :: width=5.0d-3
character(len=3) :: numb
character(len=19) :: forma
integer :: i, iwann, ipt, iipt, ityp, irpts
real(gp) :: emax, emin, inc, ener, prefac
real(gp), dimension(ndos) :: gauss
real(gp), dimension(ntypes+1,ndos) :: dos        
real(gp), dimension(ndos) :: Idos

! Calculating spread of the DOS
emax = maxval(diag)
emin = minval(diag)
inc = (emax-emin) / real(ndos-100,kind=8)  !the minus is because we want to span a bigger area than just the sprd of the energies (gaussian width)
prefac =  1.0_dp/(width*sqrt(2.0_dp*3.14159265_dp))
dos =0.0_dp

do irpts = 1, nrpts
   do iwann = 1, nwann
      ! Calculate the normalized gaussian
      do ipt = 1, ndos
         ener = emin + real(ipt-51,kind=8)*inc
         gauss(ipt) = prefac * exp(-(ener - diag(irpts,iwann))**2 / (2.0_dp*width**2))            
         ! Add the Gaussian to the DOS
         dos(1,ipt) = dos(1,ipt) + gauss(ipt)
         ityp = types(irpts,iwann) + 1
         dos(ityp,ipt) = dos(ityp,ipt) + gauss(ipt)
      end do
   end do
end do
write(numb,'(I3)') ntypes + 1
forma = '(E14.6,2x,'//numb//'E14.6)'

open(22, file=trim(filename), status='unknown')
write(22,'(A)') '# Wannier DOS file'
write(22,'(A,2x,i4)') '# Number of Wannier functions:',nwann
write(22,'(A,2x,i4)') '# Number of real space points:',nrpts
write(22,'(A,2x,i4)') '# Number of energy points:',ndos
write(22,'(A,2x,E14.6)') '# Width of the Gaussians:', width
do ipt = 1, ndos
   ener = emin + real(ipt-51,kind=8)*inc
   write(22,forma) ener, (dos(i,ipt), i=1,ntypes+1)
end do
close(22)

open(22, file='integrated_'//trim(filename), status='unknown')
write(22,'(A)') '# Wannier DOS file'
write(22,'(A,2x,i4)') '# Number of Wannier functions:',nwann
write(22,'(A,2x,i4)') '# Number of real space points:',nrpts
write(22,'(A,2x,i4)') '# Number of energy points:',ndos
write(22,'(A,2x,E14.6)') '# Width of the Gaussians:', width
Idos = 0.0_dp
do ipt = 1, ndos
   ener = emin + real(ipt-51,kind=8)*inc
   do ityp = 2, ntypes
      do iipt = 1, ipt
         Idos(ipt) = Idos(ipt) + dos(ityp,iipt)*inc
      end do
   end do
   write(22,'(E14.6,2x,E14.6)') ener, Idos(ipt)
end do
close(22)


end subroutine wannier_dos

subroutine wannier_projected_dos(filename,nrpts,nwann,norb,umn,ntypes,types,eigen)
use module_types
implicit none
character(len=*), intent(in) :: filename
integer, intent(in) :: nrpts                                               ! Number of unit cells (r-points)
integer, intent(in) :: nwann                                               ! Number of Wannier functions 
integer, intent(in) :: ntypes                                              ! Number of types of Wannier functions
integer, intent(in) :: norb                                                ! Number of orbitals
integer, dimension(nrpts,nwann), intent(in) :: types                       ! Types of the Wannier functions
real(gp), dimension(nwann,norb), intent(in) :: umn                          ! Wannier transformation matrix
real(gp), dimension(norb),intent(in) :: eigen                              ! Diagonal elements of the Hamiltonian in Wannier basis
!Local variables
integer, parameter :: ndos = 1000
real(gp), parameter :: width=1.0d-2
character(len=3) :: numb
character(len=19) :: forma
integer :: i, iwann, ipt, iipt, ityp, irpts,iorb
real(gp) :: emax, emin, inc, ener, prefac, fac
real(gp), dimension(ndos) :: gauss
!real(gp), dimension(ntypes+1,ndos) :: dos        
real(gp), dimension(nwann+1,ndos) :: dos        
real(gp), dimension(ndos) :: Idos

! Calculating spread of the DOS
emax = maxval(eigen)
emin = minval(eigen)
inc = (emax-emin) / real(ndos-100,kind=8)  !the minus is because we want to span a bigger area than just the sprd of the energies (gaussian width)
fac =  1.0_dp/(width*sqrt(2.0_dp*3.14159265_dp))
dos =0.0_dp

do irpts = 1, nrpts
   do iwann = 1, nwann
      do iorb = 1, norb
         ! Calculate the normalized gaussian
         prefac = fac * umn(iwann,iorb)**2
         do ipt = 1, ndos
            ener = emin + real(ipt-51,kind=8)*inc
            gauss(ipt) = prefac * exp(-(ener - eigen(iorb))**2 / (2.0_dp*width**2))            
            ! Add the Gaussian to the DOS
            dos(1,ipt) = dos(1,ipt) + gauss(ipt)
            ityp = iwann+1 !types(irpts,iwann) + 1
            dos(ityp,ipt) = dos(ityp,ipt) + gauss(ipt)
         end do
      end do
   end do
end do
write(numb,'(I3)') nwann +1
forma = '(E14.6,2x,'//numb//'E14.6)'

open(22, file=trim(filename), status='unknown')
write(22,'(A)') '# Wannier DOS file'
write(22,'(A,2x,i4)') '# Number of Wannier functions:',nwann
write(22,'(A,2x,i4)') '# Number of real space points:',nrpts
write(22,'(A,2x,i4)') '# Number of energy points:',ndos
write(22,'(A,2x,E14.6)') '# Width of the Gaussians:', width
do ipt = 1, ndos
   ener = emin + real(ipt-51,kind=8)*inc
!   write(22,forma) ener, (dos(i,ipt), i=1,ntypes+1)
   write(22,forma) ener, (dos(i,ipt), i=1,nwann+1)
end do
close(22)

open(22, file='integrated_'//trim(filename), status='unknown')
write(22,'(A)') '# Wannier DOS file'
write(22,'(A,2x,i4)') '# Number of Wannier functions:',nwann
write(22,'(A,2x,i4)') '# Number of real space points:',nrpts
write(22,'(A,2x,i4)') '# Number of energy points:',ndos
write(22,'(A,2x,E14.6)') '# Width of the Gaussians:', width
Idos = 0.0_dp
do ipt = 1, ndos
   ener = emin + real(ipt-51,kind=8)*inc
   do ityp = 2, nwann+1
      do iipt = 1, ipt
         Idos(ipt) = Idos(ipt) + dos(ityp,iipt)*inc
      end do
   end do
   write(22,'(E14.6,2x,E14.6)') ener, Idos(ipt)
end do
close(22)


end subroutine wannier_projected_dos

subroutine read_eigenvalues(filename,norb,nkpt,eigen)
use module_types
implicit none
character(len=*), intent(in) :: filename
integer, intent(in) :: norb, nkpt
real(gp),dimension(nkpt,norb),intent(out) :: eigen
!Local variables
logical :: file_exist
integer :: ikpt,iorb,iiorb,iikpt,ierr

inquire(file=filename,exist=file_exist)
if (.not. file_exist) then
   write(*,'(A)') 'ERROR : Input file,',filename,',not found !'
   write(*,'(A)') 'CORRECTION: Create or give correct input.inter file.'
   call mpi_finalize(ierr)
   stop
end if

open(22,file=filename,status='old')
do ikpt = 1, nkpt
   do iorb = 1, norb
      read(22,*) iiorb, iikpt, eigen(ikpt,iorb)
   end do
end do
close(22)

end subroutine read_eigenvalues

subroutine read_amn_header(filename,nproj,nband,nkpt)
use module_types
use module_interfaces
implicit none
character(len=*),intent(in) :: filename
integer, intent(out) :: nproj,nband,nkpt
!Local variables
logical :: file_exist
integer :: ierr

inquire(file=trim(filename)//'.amn',exist=file_exist)
if (.not. file_exist) then
   write(*,'(A)') 'ERROR : Input file,',trim(filename)//'.amn',', not found !'
   write(*,'(A)') 'CORRECTION: Create or give correct input.inter file.'
   call mpi_finalize(ierr)
   stop
end if

open(22,file=trim(filename)//'.amn',status='old')
read(22,*) ! skip first line which is a comment
read(22,*) nband, nkpt, nproj
close(22)
end subroutine read_amn_header

subroutine read_amn(filename,amn,nproj,nband,nkpt)
use module_types
use module_interfaces
implicit none
character(len=*),intent(in) :: filename
integer, intent(in) :: nproj, nband, nkpt
real(gp),dimension(nband,nproj), intent(out) :: amn
!Local variables
logical :: file_exist
integer :: int1, int2, int3, int4, nk, np, nb
integer :: ierr
real(gp) :: r1

inquire(file=trim(filename)//'.amn',exist=file_exist)
if (.not. file_exist) then
   write(*,'(A)') 'ERROR : Input file,',trim(filename)//'.amn',', not found !'
   write(*,'(A)') 'CORRECTION: Create or give correct input.inter file.'
   call mpi_finalize(ierr)
   stop
end if

open(22,file=trim(filename)//'.amn',status='old')
read(22,*) ! skip first line which is a comment
read(22,*) !skip this line: nband, nkpt, nproj
do nk=1, nkpt 
   do np=1, nproj
      do nb=1, nband
         read(22,*) int2, int3, int4, amn(nb,np), r1
      end do
   end do
end do
close(22)

end subroutine read_amn

!>  This routine reads an .nnkp file and returns the types of projectors
subroutine read_proj(seedname, n_kpts, n_proj, l, mr)
   implicit none
   ! I/O variables
   character(len=16),intent(in) :: seedname
   integer, intent(in) :: n_kpts, n_proj
   integer, dimension(n_proj), intent(out) :: l, mr
   ! Local variables
   integer :: i, j
   character *16 :: char1, char2, char3, char4
   logical :: calc_only_A
   real :: real_latt(3,3), recip_latt(3,3)
   real, dimension(n_kpts,3) :: kpts
   real(kind=8), dimension(n_proj,3) :: ctr_proj, x_proj, z_proj
   integer, dimension(n_proj) :: rvalue
   real, dimension(n_proj) :: zona


   OPEN(11, FILE=trim(seedname)//'.nnkp', STATUS='OLD')
   READ(11,*) ! skip first line


   !=====calc_only_A=====!
   READ(11,*) char1, char2, char3
   if (char3 .eq. 'T') then 
      calc_only_A=.TRUE.
   else 
      if (char3 .eq. 'F') then 
         calc_only_A=.FALSE.
      end if
   end if

   !=====real_lattice=====!
   READ(11,*) char1, char2 ! skip "begin real_lattice"
   READ(11,*) ((real_latt(i,j), j=1,3), i=1,3)
   READ(11,*) char3, char4 ! skip "end real_lattice"

   !=====recip_lattice=====!
   READ(11,*) char1, char2 ! skip "begin recip_lattice"
   READ(11,*) ((recip_latt(i,j), j=1,3), i=1,3)
   READ(11,*) char3, char4 ! skip "end recip_lattice"

   !=====kpoints=====!
   READ(11,*) char1, char2 ! skip "begin kpoints"
   READ(11,*) char3
   if (char3 .ne. 'end') then ! verify that there are kpoints
      BACKSPACE 11
      READ(11,*) ! skip n_kpts
      READ(11,*) ((kpts(i,j), j=1,3), i=1,n_kpts)
      READ(11,*) char3, char4 ! skip "end kpoints"
   end if

   !=====projections=====!
   READ(11,*) char1, char2 ! skip "begin projections"
   READ(11,*) char3
   if (char3 .ne. 'end') then ! verify that there are projections
      BACKSPACE 11
      READ(11,*) ! skip n_proj
      READ(11,*) ((ctr_proj(i,j), j=1,3), l(i), mr(i), rvalue(i), (z_proj(i,j), j=1,3), (x_proj(i,j), j=1,3), zona(i), i=1,n_proj)
      READ(11,*) char3, char4 ! skip "end projections"
   end if

   close(11)

END SUBROUTINE read_proj

subroutine character_list(nwann,nproj,tmatrix,plotwann,ncenters,wann_list,l,mr)
   use BigDFT_API
   use module_types
   use module_interfaces
   implicit none
   ! I/O variables
   integer, intent(in) :: nwann, nproj,plotwann
   integer, dimension(nproj), intent(in) :: l, mr
   real(gp), dimension(nwann,nproj), intent(in) :: tmatrix
   integer, dimension(plotwann), intent(in) :: ncenters
   integer, dimension(nwann),intent(in) :: wann_list
   !Local variables
   character(len=*),parameter :: subname='character_list'
   character(len=4) :: num
   character(len=17):: forma
   character(len=10), dimension(nproj) :: label
   integer :: np, np2, iwann,iiwann, ntype, ii, i_stat
   real(gp), dimension(:,:), allocatable :: Wpweight
   real(gp), dimension(:), allocatable :: norm 
   character(len=10),dimension(:), allocatable :: Wplabel
   integer, dimension(nproj) :: l_used, mr_used

   !Start by finding the projector labels
   do np = 1, nproj
      if (l(np)==0) then   ! s orbital
         label(np) = 's'
      end if
      if (l(np)==1) then   ! p orbitals
         if (mr(np)==1) label(np) = 'pz'
         if (mr(np)==2) label(np) = 'px'
         if (mr(np)==3) label(np) = 'py'
      end if
      if (l(np)==2) then   ! d orbitals
         if (mr(np)==1) label(np) = 'dz2'
         if (mr(np)==2) label(np) = 'dxz'
         if (mr(np)==3) label(np) = 'dyz'
         if (mr(np)==4) label(np) = 'dx2-y2'
         if (mr(np)==5) label(np) = 'dxy'
      endif
      if (l(np)==3) then   ! f orbitals
         if (mr(np)==1) label(np) = 'fz3'  
         if (mr(np)==2) label(np) = 'fxz2'
         if (mr(np)==3) label(np) = 'fyz2'
         if (mr(np)==4) label(np) = 'fz(x2-y2)'
         if (mr(np)==5) label(np) = 'fxyz'
         if (mr(np)==6) label(np) = 'fx(x2-3y2)'
         if (mr(np)==7) label(np) = 'fy(3x2-y2)'
      endif
      if (l(np)==-1) then  !  sp hybrids
         if (mr(np)==1) label(np) = 'sp-1' 
         if (mr(np)==2) label(np) = 'sp-2'
      end if
      if (l(np)==-2) then  !  sp2 hybrids 
         if (mr(np)==1) label(np) = 'sp2-1' 
         if (mr(np)==2) label(np) = 'sp2-2'
         if (mr(np)==3) label(np) = 'sp2-3'
      end if
      if (l(np)==-3) then  !  sp3 hybrids
         if (mr(np)==1) label(np) = 'sp3-1'
         if (mr(np)==2) label(np) = 'sp3-2'
         if (mr(np)==3) label(np) = 'sp3-3'
         if (mr(np)==4) label(np) = 'sp3-4'
      end if
      if (l(np)==-4) then  !  sp3d hybrids
         if (mr(np)==1) label(np) = 'sp3d-1'
         if (mr(np)==2) label(np) = 'sp3d-2'
         if (mr(np)==3) label(np) = 'sp3d-3'
         if (mr(np)==4) label(np) = 'sp3d-4'
         if (mr(np)==5) label(np) = 'sp3d-5'
      end if
      if (l(np)==-5) then  ! sp3d2 hybrids
         if (mr(np)==1) label(np) = 'sp3d2-1'
         if (mr(np)==2) label(np) = 'sp3d2-2'
         if (mr(np)==3) label(np) = 'sp3d2-3'
         if (mr(np)==4) label(np) = 'sp3d2-4'
         if (mr(np)==5) label(np) = 'sp3d2-5'
         if (mr(np)==6) label(np) = 'sp3d2-6'
      end if
   end do   

   ! count the number of projector types
   ntype = 0
   l_used = -99
   mr_used =-99
   loop_np1: do np = 1, nproj
      !Check wether we have already used this projector type
      do np2 = 1,nproj
         if(l(np) == l_used(np2) .and. mr(np) == mr_used(np2)) then
            cycle loop_np1
         end if
      end do
      ntype=ntype+1
      l_used(ntype) = l(np)
      mr_used(ntype) = mr(np)
   end do loop_np1

   allocate(Wpweight(nwann,ntype),stat=i_stat)
   call memocc(i_stat,Wpweight,'Wpweight',subname)
   allocate(Wplabel(ntype),stat=i_stat)
   call memocc(i_stat,Wplabel,'Wplabel',subname)

   ! Construct the weights of each type
   ii = 0
   l_used = -99
   mr_used =-99
   Wpweight = 0.0d0
   Wplabel = 'unspec'
   loop_np: do np = 1, nproj
      !Check wether we have already used this projector type
      do np2 = 1,nproj
         if(l(np) == l_used(np2) .and. mr(np) == mr_used(np2)) then
            cycle loop_np
         end if
      end do
      ii=ii+1
      l_used(ii) = l(np)
      mr_used(ii) = mr(np)
      do np2 = 1, nproj
         if(l(np2) == l(np) .and. mr(np2) == mr(np)) then
            do iwann = 1, nwann
              Wpweight(iwann,ii) = Wpweight(iwann,ii) + tmatrix(iwann,np2)**2
              Wplabel(ii) = label(np)
            end do
         end if
      end do
   end do loop_np

   allocate(norm(nwann), stat=i_stat)
   call memocc(i_stat,norm,'norm',subname)

   !calcualte norm
   norm = 0.0d0
   do iwann=1,nwann
      do np=1,ntype
         norm(iwann) = norm(iwann) + Wpweight(iwann,np)**2
      end do
      norm(iwann) = sqrt(norm(iwann))
   end do

   ! Print the information
!   if(iproc==0) then
     write(*,*) 'Analysis of the symmetry types of the Wannier functions'
     write(num,'(I4)') ntype
     forma = '(23x,'//trim(num)//'(A,17x))'
     write(*,trim(forma))(Wplabel(ii),ii=1,ntype)
     do iwann = 1, plotwann
        iiwann = wann_list(iwann)
        write(*,*) iiwann, iwann, (Wpweight(iiwann,ii)/norm(iiwann), ii=1,ntype)
     end do
!   end if

end subroutine character_list

!!$subroutine build_stereographic_graph(natoms,proj,nsurf,ncenters,Zatoms,normal,NeglectPoint)
!!$   use module_interfaces
!!$   use module_types
!!$   implicit none
!!$   integer, intent(in) :: natoms, nsurf, NeglectPoint
!!$   real(gp), dimension(natoms,3), intent(in) :: proj
!!$   integer, dimension(nsurf), intent(in) :: ncenters
!!$   integer, dimension(natoms,nsurf),intent(in) :: Zatoms
!!$   real(gp), dimension(3), intent(in) :: normal
!!$   ! Local variables
!!$   integer :: isurf,ipts,i, nsurftot,npts,ndecimal,ncent
!!$   character(len=14) :: surfname
!!$   character(len=20) :: forma
!!$   logical :: condition
!!$
!!$
!!$   open(22,file='proj.surf', status='unknown')
!!$
!!$   !First line is just a comment 
!!$   write(22,'(A)') 'Stereographic graph'
!!$
!!$   !Write the dimensions of the box (simple cubic)
!!$   write(22,'(E14.6, 2x, E14.6, 2x, E14.6)') maxval(proj(:,1))-minval(proj(:,1)), 0.0, maxval(proj(:,2))-minval(proj(:,2))  !dxx dyx dyy
!!$   write(22,'(E14.6, 2x, E14.6, 2x, E14.6)') 0.0, 0.0,  maxval(proj(:,3))-minval(proj(:,3))                                   !dzx dzy dzz
!!$
!!$   nsurftot = 0
!!$   npts = 0
!!$   do isurf = 1, nsurf
!!$      !MUST eliminate the wanniers that are less then 3 centers
!!$      if(ncenters(isurf) < 3) cycle
!!$      ! or the 3 centers containing the neglected point
!!$      condition = (Zatoms(1,isurf) == NeglectPoint .or. Zatoms(2,isurf) == NeglectPoint .or. Zatoms(3,isurf) == NeglectPoint)
!!$      if(ncenters(isurf) == 3 .and. condition) cycle
!!$      nsurftot = nsurftot + 1
!!$      ! Must eliminate the neglected point from other surfaces
!!$      condition = .false.
!!$      do ipts = 1, ncenters(isurf)
!!$         condition = condition .or. Zatoms(ipts,isurf) == NeglectPoint
!!$      end do
!!$      if(condition) then
!!$         npts = npts + ncenters(isurf)-1
!!$      else
!!$         npts = npts + ncenters(isurf)
!!$      end if
!!$   end do
!!$   !Fourth line: number of surfaces, total num_polys, total num_points
!!$   write(22,'(I4, 2x, I4, 2x, I4)') nsurftot, nsurftot, npts
!!$
!!$   do isurf = 1, nsurf
!!$      !MUST eliminate the wanniers that are less then 3 centers
!!$      if(ncenters(isurf) < 3) cycle
!!$      ! or the 3 centers containing the neglected point
!!$      condition = (Zatoms(1,isurf) == NeglectPoint .or. Zatoms(2,isurf) == NeglectPoint .or. Zatoms(3,isurf) == NeglectPoint)
!!$      if(ncenters(isurf) == 3 .and. condition) cycle
!!$
!!$      ! Determining the format (number of integers) for surface name 
!!$      ndecimal = 1
!!$      ncent = ncenters(isurf)
!!$      do
!!$        if(int(ncent / 10) == 0) exit
!!$        ncent = ncent / 10
!!$        ndecimal = ndecimal + 1
!!$      end do
!!$      write(forma,'(I1)') ndecimal
!!$      forma = '(I'//trim(forma)//')'
!!$
!!$      !Name of the surface
!!$      write(surfname,forma) ncenters(isurf)
!!$      surfname = '2e - '//trim(surfname)//'centers'
!!$      write(22,'(A)') trim(surfname)
!!$
!!$      !num_polys and num_points (Must eliminate the neglected point from other surfaces)
!!$      condition = .false.
!!$      do ipts = 1, ncenters(isurf)
!!$         condition = condition .or. Zatoms(ipts,isurf) == NeglectPoint
!!$      end do
!!$      if(condition) then
!!$         write(22,'(I4, 2x, I4)') 1, ncenters(isurf)-1
!!$         !number of vertices, i_1 i_2 i_3 ... i_n (index of the vertices)
!!$         write(forma,'(I4)') ncenters(isurf) - 1
!!$         forma = '(I4,'//trim(forma)//'(2x,I4))'
!!$         write(22,trim(forma)) ncenters(isurf) - 1, (i,i=1,ncenters(isurf) - 1 )!(Zatoms(i,isurf),i=1,ncenters(isurf))
!!$      else
!!$         write(22,'(I4, 2x, I4)') 1, ncenters(isurf)
!!$         !number of vertices, i_1 i_2 i_3 ... i_n (index of the vertices)
!!$         write(forma,'(I4)') ncenters(isurf)
!!$         forma = '(I4,'//trim(forma)//'(2x,I4))'
!!$         write(22,trim(forma)) ncenters(isurf), (i,i=1,ncenters(isurf))!(Zatoms(i,isurf),i=1,ncenters(isurf))
!!$      end if
!!$
!!$      do ipts = 1, ncenters(isurf)
!!$         if(Zatoms(ipts,isurf) == NeglectPoint) cycle
!!$         !coordinates of the vertices (x y z) and the normal to the surface at these vertices (nx ny nz)
!!$         write(22,'(5(E14.6, 2x),E14.6)') proj(Zatoms(ipts,isurf),1), proj(Zatoms(ipts,isurf),2), proj(Zatoms(ipts,isurf),3),&
!!$              (normal(i), i=1,3)
!!$      end do
!!$   end do
!!$
!!$end subroutine build_stereographic_graph


!!$subroutine output_stereographic_graph(natoms,proj,nsurf,ncenters,Zatoms,nvertex,vertex,npoly,poly,normal,NeglectPoint) 
!!$   use module_interfaces
!!$   use module_types
!!$   implicit none
!!$   integer, intent(in) :: natoms, nsurf, NeglectPoint
!!$   real(gp), dimension(natoms,3), intent(in) :: proj
!!$   integer, dimension(nsurf), intent(in) :: ncenters, npoly, nvertex
!!$   integer, dimension(natoms,nsurf) :: Zatoms                          ! indexes of all the atoms spanned by the Wannier function
!!$   integer, dimension(nsurf,maxval(npoly),3),intent(in) :: poly
!!$   integer, dimension(maxval(nvertex),nsurf), intent(in) :: vertex     ! indexes of the vertex (on the convex hull) of the surface
!!$   real(gp), dimension(3), intent(in) :: normal
!!$   ! Local variables
!!$   integer :: isurf,ipts,i, nsurftot,npts,ndecimal,ncent, num_poly
!!$   character(len=14) :: surfname   
!!$   character(len=20) :: forma
!!$   logical :: condition
!!$   
!!$   
!!$   open(22,file='proj.surf', status='unknown')
!!$   
!!$   !First line is just a comment 
!!$   write(22,'(A)') 'Stereographic graph'
!!$   
!!$   !Write the dimensions of the box (simple cubic)
!!$   write(22,'(E14.6, 2x, E14.6, 2x, E14.6)') maxval(proj(:,1))-minval(proj(:,1)), 0.0, maxval(proj(:,2))-minval(proj(:,2))  !dxx dyx dyy
!!$   write(22,'(E14.6, 2x, E14.6, 2x, E14.6)') 0.0, 0.0,  maxval(proj(:,3))-minval(proj(:,3))                                   !dzx dzy dzz
!!$
!!$   nsurftot = 0
!!$   npts = 0
!!$   num_poly = 0
!!$   do isurf = 1, nsurf
!!$      !MUST eliminate the wanniers that are less then 3 centers
!!$      if(ncenters(isurf) < 3) cycle
!!$      ! or the 3 centers containing the neglected point
!!$      condition = (Zatoms(1,isurf) == NeglectPoint .or. Zatoms(2,isurf) == NeglectPoint .or. Zatoms(3,isurf) == NeglectPoint)
!!$      if(ncenters(isurf) == 3 .and. condition) cycle
!!$      nsurftot = nsurftot + 1
!!$      ! Must eliminate the neglected point from other surfaces
!!$      do ipts = 1, npoly(isurf)
!!$         condition = .false.
!!$         do i = 1, 3
!!$            condition = condition .or. poly(ipts,i,isurf) == NeglectPoint
!!$         end do
!!$         if(.not. condition) num_poly = num_poly + 1
!!$      end do
!!$      if(condition) then
!!$         npts = npts + nvertex(isurf)-1
!!$      else
!!$         npts = npts + nvertex(isurf)
!!$      end if
!!$   end do
!!$
!!$   !Fourth line: number of surfaces, total num_polys, total num_points
!!$   write(22,'(I4, 2x, I4, 2x, I4)') nsurftot, num_poly, npts
!!$
!!$   do isurf = 1, nsurf
!!$      !MUST eliminate the wanniers that are less then 3 centers
!!$      if(ncenters(isurf) < 3) cycle
!!$      ! or the 3 centers containing the neglected point
!!$      condition = (Zatoms(1,isurf) == NeglectPoint .or. Zatoms(2,isurf) == NeglectPoint .or. Zatoms(3,isurf) == NeglectPoint)
!!$      if(ncenters(isurf) == 3 .and. condition) cycle
!!$
!!$      ! Determining the format (number of integers) for surface name 
!!$      ndecimal = 1
!!$      ncent = ncenters(isurf)
!!$      do
!!$        if(int(ncent / 10) == 0) exit
!!$        ncent = ncent / 10
!!$        ndecimal = ndecimal + 1
!!$      end do
!!$      write(forma,'(I1)') ndecimal
!!$      forma = '(I'//trim(forma)//')'
!!$
!!$      !Name of the surface
!!$      write(surfname,forma) ncenters(isurf)
!!$      surfname = '2e - '//trim(surfname)//'centers'
!!$      write(22,'(A)') trim(surfname)
!!$
!!$      !num_polys and num_points (Must eliminate the neglected point from other surfaces)
!!$      condition = .false.
!!$      num_poly = 0
!!$      do ipts = 1, npoly(isurf)
!!$         condition = .false.
!!$         do i = 1, 3
!!$            condition = condition .or. poly(ipts,i,isurf) == NeglectPoint
!!$         end do
!!$         if(.not. condition)num_poly = num_poly + 1
!!$      end do
!!$      if(condition) then
!!$         write(22,'(I4, 2x, I4)') num_poly, nvertex(isurf)-1
!!$         !number of vertices, i_1 i_2 i_3 ... i_n (index of the vertices)
!!$         do ipts = 1, npoly(isurf)
!!$            condition = .false.
!!$            do i = 1, 3
!!$               condition = condition .or. poly(ipts,i,isurf) == NeglectPoint
!!$            end do
!!$            if(condition) cycle
!!$            write(forma,'(I4)') nvertex(isurf) - 1
!!$            forma = '(I4,'//trim(forma)//'(2x,I4))'
!!$            write(22,trim(forma)) 3, (poly(ipts,i,isurf),i=1,3) ! Only triangles
!!$         end do
!!$      else
!!$         write(22,'(I4, 2x, I4)') num_poly, nvertex(isurf)
!!$         !number of vertices, i_1 i_2 i_3 ... i_n (index of the vertices)
!!$         do ipts = 1, npoly(isurf)
!!$            condition = .false.
!!$            do i = 1, 3
!!$               condition = condition .or. poly(ipts,i,isurf) == NeglectPoint
!!$            end do
!!$            if(condition) cycle
!!$            write(forma,'(I4)') ncenters(isurf)
!!$            forma = '(I4,'//trim(forma)//'(2x,I4))'
!!$            write(22,trim(forma)) 3, (poly(ipts,i,isurf),i=1,3) ! Only triangles
!!$         end do
!!$      end if
!!$
!!$      do ipts = 1, nvertex(isurf)
!!$         if(vertex(ipts,isurf) == NeglectPoint) cycle
!!$         !coordinates of the vertices (x y z) and the normal to the surface at these vertices (nx ny nz)
!!$         write(22,'(5(E14.6, 2x),E14.6)') proj(vertex(ipts,isurf),1), proj(vertex(ipts,isurf),2), proj(vertex(ipts,isurf),3),&
!!$              (normal(i), i=1,3)
!!$      end do
!!$   end do
!!$
!!$end subroutine output_stereographic_graph


!!$subroutine build_stereographic_surface(nsurf,atoms,rxyz,ncenters,list)
!!$   use module_types
!!$   implicit none
!!$   integer, intent(in) :: nsurf
!!$   type(atoms_data), intent(in) :: atoms
!!$   real(gp), dimension(3,atoms%nat), intent(in) :: rxyz
!!$   integer, dimension(nsurf), intent(in) :: ncenters
!!$   integer, dimension(atoms%nat,nsurf),intent(in) :: list
!!$   !Local variables
!!$   logical :: inverted
!!$   integer :: iat, i, j, isurf, ifacet, NeglectPoint, SNeglectPoint
!!$   integer :: npts, nfacets, ifcts, iisurf, ipoly, numdim
!!$   real(gp) :: rad, dist, norm
!!$   real(gp), dimension(atoms%nat,3) :: projN, projS      ! atom positions in the projections
!!$   real(gp), allocatable :: proj(:,:)
!!$   real(gp), dimension(3) :: CM, refpos, normal, southref, newpt
!!$   integer, dimension(nsurf) :: npoly
!!$   integer, dimension(2,3) :: newfct
!!$   integer, allocatable :: facets(:,:), newFacets(:,:), tmp2(:,:), tmp3 (:,:,:), poly(:,:,:)
!!$   integer, allocatable :: nvertex(:), vertex(:,:)
!!$   real(gp), allocatable :: tmp(:,:)
!!$
!!$!!$ Should change this
!!$   ! For now choosing the reference point by hand
!!$    refpos(1) = rxyz(1,78) ; refpos(2) = rxyz(2,78) ; refpos(3) = rxyz(3,78)
!!$!   refpos(1) = (rxyz(1,22) + rxyz(1,43) + rxyz(1,39) + rxyz(1,26) + rxyz(1,6)) / 5.0
!!$!   refpos(2) = (rxyz(2,22) + rxyz(2,43) + rxyz(2,39) + rxyz(2,26) + rxyz(2,6)) / 5.0
!!$!   refpos(3) = (rxyz(3,22) + rxyz(3,43) + rxyz(3,39) + rxyz(3,26) + rxyz(3,6)) / 5.0
!!$!!$ END SHOULD CHANGE THIS
!!$
!!$   open(22,file='pos_ref.xyz', status='unknown')
!!$   write(22,'(I4)') atoms%nat+1
!!$   write(22,*) !skip this line
!!$   write(22,'(A,3(2x,E14.6))'), 'X',(refpos(i),i=1,3)
!!$   do i = 1, atoms%nat
!!$      write(22,'(A,3(2x,E14.6))'),atoms%atomnames(atoms%iatype(i)),rxyz(1,i),rxyz(2,i),rxyz(3,i) 
!!$   end do
!!$   close(22)
!!$
!!$   ! Calculate Center of mass
!!$   CM = 0.0
!!$   do iat = 1, atoms%nat
!!$      do j = 1, 3
!!$         CM(j) = CM(j) + rxyz(j,iat) / atoms%nat
!!$      end do
!!$   end do
!!$
!!$   !Calculate the radius of the sphere (choose it to be the biggest distance from the CM)
!!$   rad= 0.0
!!$   do iat = 1, atoms%nat
!!$      dist = 0.0
!!$      do j = 1, 3
!!$         dist = dist + (rxyz(j,iat) - CM(j))**2
!!$      end do
!!$      rad = max(dist, rad)
!!$   end do
!!$   rad =sqrt(rad)
!!$
!!$
!!$   ! calculate northern projection
!!$   call stereographic_projection(atoms%nat, rxyz, refpos, CM, rad, projN, normal, NeglectPoint)
!!$   
!!$   ! Copy projN to proj
!!$   allocate(proj(atoms%nat, 3))
!!$   do iat = 1, atoms%nat
!!$      do i = 1, 3
!!$         proj(iat,i) = projN(iat,i)
!!$      end do
!!$   end do
!!$
!!$   call write_stereographic_projection(22, 'proj.xyz    ', atoms, proj, NeglectPoint)
!!$ 
!!$   ! Calculate the southern pole
!!$   do i = 1, 3
!!$      southref(i) = (CM(i) - refpos(i))
!!$   end do
!!$   norm = sqrt(southref(1)**2 + southref(2)**2 + southref(3)**2 )
!!$   do i = 1, 3
!!$      southref(i) = southref(i)*rad/norm + CM(i)
!!$   end do
!!$
!!$   ! Calculate the southern projection
!!$   call stereographic_projection(atoms%nat, rxyz, southref, CM, rad, projS, normal, SNeglectPoint)
!!$
!!$
!!$   open(22,file='Spos_ref.xyz', status='unknown')
!!$   write(22,'(I4)') atoms%nat+1
!!$   write(22,*) !skip this line
!!$   write(22,'(A,3(2x,E14.6))'), 'X',(southref(i),i=1,3)
!!$   do i = 1, atoms%nat
!!$      write(22,'(A,3(2x,E14.6))'),atoms%atomnames(atoms%iatype(i)),rxyz(1,i),rxyz(2,i),rxyz(3,i) 
!!$   end do
!!$   close(22)
!!$
!!$   call write_stereographic_projection(22, 'Sproj.xyz   ', atoms, projS, SNeglectPoint)
!!$
!!$   allocate(nvertex(nsurf)) 
!!$   allocate(vertex(maxval(ncenters),nsurf))
!!$   vertex = 0
!!$
!!$   npts = atoms%nat
!!$   do isurf = 1, nsurf
!!$
!!$      numdim = ncenters(isurf)*(ncenters(isurf)-1)*(ncenters(isurf)-2)/6
!!$      if(allocated(facets))deallocate(facets)
!!$      allocate(facets(numdim,3))
!!$ 
!!$!      call convex_hull_construction_3D_CS1989(natoms,rxyz,ncenters,list,nvertex,vertex,nfacets,facets)
!!$      ! Don't really need the convex hull (only worth it to eliminate points in the interior)
!!$      ! Make all the possible triangles, without permutations (n!/3!(n-3)!)
!!$      call make_facets(ncenters(isurf),list(1,isurf),nvertex(isurf),vertex(1,isurf),nfacets,facets)
!!$
!!$      ! Copy facets to newFacets
!!$      if(allocated(newFacets))deallocate(newFacets)
!!$      allocate(newFacets(nfacets, 3))
!!$      do iat = 1, nfacets
!!$         do i = 1, 3
!!$            newFacets(iat,i) = facets(iat,i)
!!$         end do
!!$      end do
!!$   
!!$      ! For the triangles not completly in the northern hemisphere
!!$      ! Calculate the interior normal of the edges of the triangles for north hemisphere projection (reference point)
!!$      ! Calculate the interior normal of the edges of the triangles for south hemisphere projection
!!$      ! Compare the two, if normal does not inverts itself(change of pole should induce a 180 rotation), should divide this edge in two (recursif)
!!$      npoly(isurf) = nfacets
!!$      do ifacet = 1, nfacets
!!$         call build_correct_triangles_for_projection(atoms%nat, rxyz, projN, projS, facets(ifacet,1),&
!!$              refpos, southref, CM, rad, inverted, newpt)
!!$
!!$         newfct(1,1) = npts + 1; newfct(1,2) = facets(ifacet,1) ; newfct(1,3) = facets(ifacet,3)
!!$         newfct(2,1) = npts + 1; newfct(2,2) = facets(ifacet,2) ; newfct(2,3) = facets(ifacet,3) 
!!$   
!!$         if(inverted) then
!!$           ! Must add the extra point to the atom positions in the projection
!!$           ! Always add at the end, such that we do not change the previous indexes in facets
!!$           if(allocated(tmp)) deallocate(tmp)
!!$           allocate(tmp(npts,3))
!!$           do iat = 1, npts 
!!$              do i = 1, 3
!!$                 tmp(iat,i) = proj(iat,i)
!!$              end do
!!$           end do
!!$           if(allocated(proj)) deallocate(proj)
!!$           allocate(proj(npts+1,3))
!!$           do iat = 1, npts 
!!$              do i = 1, 3
!!$                 proj(iat,i) = tmp(iat,i)
!!$              end do
!!$           end do
!!$           do i = 1, 3
!!$              proj(npts + 1,i) = newpt(i)
!!$           end do
!!$           npts = npts + 1
!!$   
!!$           ! Must change facets accordingly
!!$           if(allocated(tmp2)) deallocate(tmp2)
!!$           allocate(tmp2(npoly(isurf),3))
!!$           do iat = 1, npoly(isurf)
!!$              do i = 1, 3
!!$                 tmp2(iat,i) = newFacets(iat,i)
!!$              end do
!!$           end do
!!$           if(allocated(newfacets)) deallocate(newfacets)
!!$           allocate(newFacets(npoly(isurf)+1,3))
!!$           i = 0
!!$           do ifcts = 1, npoly(isurf)
!!$              if(ifcts == ifacet) cycle
!!$              i = i + 1
!!$              do j = 1, 3
!!$                 newFacets(i,j) = tmp2(ifcts,j)
!!$              end do
!!$           end do
!!$           ! Add the two new facets
!!$           do i = 0 , 1
!!$              do j = 1, 3 
!!$                 newFacets(npoly(isurf)+i,j) = newfct(1+i,j)   
!!$              end do
!!$           end do
!!$           npoly(isurf) = npoly(isurf) + 1
!!$         end if
!!$      end do  ! loop on facets
!!$
!!$      ! Must conserve the information for each surface
!!$      if (isurf >= 2) then
!!$         ! To do this, check that npoly is not greater then the one used to allocate
!!$         ! If it is the case, we must resize the array
!!$         if(npoly(isurf) > maxval(npoly(1:isurf-1))) then
!!$            if(allocated(tmp3)) deallocate(tmp3)
!!$            allocate(tmp3(nsurf,maxval(npoly(1:isurf-1)),3))
!!$            tmp3 = 0.0
!!$            do iisurf = 1, isurf-1
!!$               do ipoly = 1, maxval(npoly(1:isurf-1))
!!$                  do i = 1, 3
!!$                     tmp3(iisurf,ipoly,i) = poly(iisurf,ipoly,i)
!!$                  end do
!!$               end do
!!$            end do 
!!$            if(allocated(poly)) deallocate(poly)
!!$            allocate(poly(nsurf,npoly(isurf),3))
!!$            poly = 0.0
!!$            do iisurf = 1, isurf-1
!!$               do ipoly = 1, maxval(npoly(1:isurf-1))
!!$                  do i = 1, 3
!!$                     poly(iisurf,ipoly,i) = tmp3(iisurf,ipoly,i)
!!$                  end do
!!$               end do
!!$            end do 
!!$            do ipoly = 1, npoly(isurf)
!!$               do i = 1, 3
!!$                  poly(isurf,ipoly,i) = newFacets(ipoly,i)
!!$               end do
!!$            end do
!!$         else !just add the information at the correct place
!!$            do ipoly = 1, npoly(isurf)
!!$               do i = 1, 3
!!$                  poly(isurf,ipoly,i) = newFacets(ipoly,i)
!!$               end do
!!$            end do
!!$         end if
!!$      else
!!$          allocate(poly(nsurf,npoly(isurf),3))
!!$          poly = 0.0
!!$          do ipoly = 1, npoly(isurf)
!!$             do i = 1, 3
!!$                poly(isurf,ipoly,i) = newFacets(ipoly,i)
!!$             end do
!!$          end do
!!$      end if
!!$   end do ! loop on surfaces
!!$   
!!$   ! Output the surfaces
!!$   call build_stereographic_graph(npts,proj,nsurf,ncenters,nvertex,vertex,npoly,poly,normal,NeglectPoint)
!!$
!!$end subroutine build_stereographic_surface


!!$subroutine output_stereographic_graph(natoms,proj,nsurf,ncenters,Zatoms,nvertex,vertex,npoly,poly,normal,NeglectPoint)
!!$   use module_interfaces
!!$   use module_types
!!$   implicit none
!!$   integer, intent(in) :: natoms, nsurf, NeglectPoint
!!$   real(gp), dimension(natoms,3), intent(in) :: proj
!!$   integer, dimension(nsurf), intent(in) :: ncenters, npoly, nvertex
!!$   integer, dimension(natoms,nsurf) :: Zatoms                          ! indexes of all the atoms spanned by the Wannier function
!!$   integer, dimension(nsurf,maxval(npoly),3),intent(in) :: poly
!!$   integer, dimension(maxval(nvertex),nsurf), intent(in) :: vertex     ! indexes of the vertex (on the convex hull) of the surface
!!$   real(gp), dimension(3), intent(in) :: normal
!!$   ! Local variables
!!$   integer :: isurf,ipts,i, nsurftot,npts,ndecimal,ncent, num_poly
!!$   character(len=14) :: surfname   
!!$   character(len=20) :: forma
!!$   logical :: condition
!!$   
!!$   
!!$   open(22,file='proj.surf', status='unknown')
!!$   
!!$   !First line is just a comment 
!!$   write(22,'(A)') 'Stereographic graph'
!!$   
!!$   !Write the dimensions of the box (simple cubic)
!!$   write(22,'(E14.6, 2x, E14.6, 2x, E14.6)') maxval(proj(:,1))-minval(proj(:,1)), 0.0, maxval(proj(:,2))-minval(proj(:,2))  !dxx dyx dyy
!!$   write(22,'(E14.6, 2x, E14.6, 2x, E14.6)') 0.0, 0.0,  maxval(proj(:,3))-minval(proj(:,3))                                   !dzx dzy dzz
!!$   
!!$   nsurftot = 0
!!$   npts = 0
!!$   num_poly = 0
!!$   do isurf = 1, nsurf
!!$      !MUST eliminate the wanniers that are less then 3 centers
!!$      if(ncenters(isurf) < 3) cycle
!!$      ! or the 3 centers containing the neglected point
!!$      condition = (Zatoms(1,isurf) == NeglectPoint .or. Zatoms(2,isurf) == NeglectPoint .or. Zatoms(3,isurf) == NeglectPoint)
!!$      if(ncenters(isurf) == 3 .and. condition) cycle
!!$      nsurftot = nsurftot + 1
!!$      ! Must eliminate the neglected point from other surfaces
!!$      do ipts = 1, npoly(isurf)
!!$         condition = .false.
!!$         do i = 1, 3
!!$            condition = condition .or. poly(ipts,i,isurf) == NeglectPoint
!!$         end do
!!$         if(.not. condition) num_poly = num_poly + 1 
!!$      end do
!!$      if(condition) then
!!$         npts = npts + nvertex(isurf)-1
!!$      else
!!$         npts = npts + nvertex(isurf)
!!$      end if
!!$   end do
!!$
!!$   !Fourth line: number of surfaces, total num_polys, total num_points
!!$   write(22,'(I4, 2x, I4, 2x, I4)') nsurftot, num_poly, npts
!!$   
!!$   do isurf = 1, nsurf
!!$      !MUST eliminate the wanniers that are less then 3 centers
!!$      if(ncenters(isurf) < 3) cycle
!!$      ! or the 3 centers containing the neglected point
!!$      condition = (Zatoms(1,isurf) == NeglectPoint .or. Zatoms(2,isurf) == NeglectPoint .or. Zatoms(3,isurf) == NeglectPoint)
!!$      if(ncenters(isurf) == 3 .and. condition) cycle
!!$
!!$      ! Determining the format (number of integers) for surface name 
!!$      ndecimal = 1
!!$      ncent = ncenters(isurf)
!!$      do
!!$        if(int(ncent / 10) == 0) exit
!!$        ncent = ncent / 10
!!$        ndecimal = ndecimal + 1
!!$      end do
!!$      write(forma,'(I1)') ndecimal
!!$      forma = '(I'//trim(forma)//')'
!!$
!!$      !Name of the surface
!!$      write(surfname,forma) ncenters(isurf)
!!$      surfname = '2e - '//trim(surfname)//'centers'
!!$      write(22,'(A)') trim(surfname)
!!$   
!!$      !num_polys and num_points (Must eliminate the neglected point from other surfaces)
!!$      condition = .false.
!!$      num_poly = 0
!!$      do ipts = 1, npoly(isurf)
!!$         condition = .false.
!!$         do i = 1, 3
!!$            condition = condition .or. poly(ipts,i,isurf) == NeglectPoint
!!$         end do
!!$         if(.not. condition)num_poly = num_poly + 1
!!$      end do
!!$      if(condition) then
!!$         write(22,'(I4, 2x, I4)') num_poly, nvertex(isurf)-1
!!$         !number of vertices, i_1 i_2 i_3 ... i_n (index of the vertices)
!!$         do ipts = 1, npoly(isurf)
!!$            condition = .false.
!!$            do i = 1, 3
!!$               condition = condition .or. poly(ipts,i,isurf) == NeglectPoint
!!$            end do
!!$            if(condition) cycle
!!$            write(forma,'(I4)') nvertex(isurf) - 1
!!$            forma = '(I4,'//trim(forma)//'(2x,I4))'
!!$            write(22,trim(forma)) 3, (poly(ipts,i,isurf),i=1,3) ! Only triangles
!!$         end do
!!$      else
!!$         write(22,'(I4, 2x, I4)') num_poly, nvertex(isurf)
!!$         !number of vertices, i_1 i_2 i_3 ... i_n (index of the vertices)
!!$         do ipts = 1, npoly(isurf)
!!$            condition = .false.
!!$            do i = 1, 3
!!$               condition = condition .or. poly(ipts,i,isurf) == NeglectPoint
!!$            end do
!!$            if(condition) cycle
!!$            write(forma,'(I4)') ncenters(isurf)
!!$            forma = '(I4,'//trim(forma)//'(2x,I4))'
!!$            write(22,trim(forma)) 3, (poly(ipts,i,isurf),i=1,3) ! Only triangles
!!$         end do
!!$      end if
!!$
!!$      do ipts = 1, nvertex(isurf)
!!$         if(vertex(ipts,isurf) == NeglectPoint) cycle
!!$         !coordinates of the vertices (x y z) and the normal to the surface at these vertices (nx ny nz)
!!$         write(22,'(5(E14.6, 2x),E14.6)') proj(vertex(ipts,isurf),1), proj(vertex(ipts,isurf),2), proj(vertex(ipts,isurf),3),&
!!$              (normal(i), i=1,3)
!!$      end do
!!$   end do
!!$
!!$end subroutine output_stereographic_graph

!!$
!!$
!!$subroutine convex_hull_construction_3D_CS1989(natoms,rxyz,ncenters,list,nvertex,vertex,nfacets,facets)
!!$   use module_interfaces
!!$   use module_types
!!$   implicit none
!!$   integer, intent(in) :: natoms, ncenters
!!$   integer, dimension(ncenters),intent(in) :: list
!!$   real(gp), dimension(3,natoms), intent(in) :: rxyz
!!$   integer, intent(out) :: nvertex,nfacets
!!$   integer, dimension(ncenters) :: vertex
!!$   integer, dimension(999,3), intent(out) :: facets
!!$   ! Local variables
!!$   logical :: allcoplanar, collinear, allcollinear
!!$   integer :: i, j, stat,
!!$
!!$
!!$   !Only build the correct facet for ncenters == 3
!!$   if(ncenters < 2) then
!!$      nfacets = 0
!!$      return
!!$   else if(ncenters == 3) then
!!$      nfacets = 1
!!$      do i = 1, 3
!!$         facets(1,i) = list(i)
!!$      end do 
!!$      return
!!$   end if
!!$
!!$   !1) Construct a tetrahedron by connecting first four points
!!$   !a) Initialize inipts
!!$   do i = 1, 4
!!$      do j = 1, 3
!!$         inipts(j,i) = rxyz(j,list(i))
!!$      end do
!!$    end do
!!$
!!$   !b) Test coplanarity
!!$   !   To do this check that volume of tetrahedron in not zero
!!$   call volume_tetrahedron(inipts,volume,stat)
!!$
!!$   if(stat == 0) 
!!$   else if(stat == 1) then
!!$      ! Coplanar because last point is inside the plane spanned by first three
!!$      if(ncenters > 4) then ! else we search for a tetrahedron by going through the list of points
!!$         do i = 5, ncenters
!!$            do j = 1, 3
!!$               initpts(4,j) = rxyz(j,list(i))
!!$            end do
!!$            call volume_tetrahedron(inipts,volume,stat)
!!$            if(stat == 0) exit  !we have found our point
!!$            if(stat .ne. 0 and i == ncenters) allcoplanar = .true.
!!$         end do
!!$      end if
!!$      if(ncenters == 4 .or. allcoplanar) then
!!$         ! if there is only four points, use graham algorithm 2D
!!$         call graham_convex_hull()
!!$      end if
!!$   else if(stat == 2) then
!!$      ! Coplanar because first three points are collinear
!!$      if(ncenters > 4) then !change the third point
!!$         do i = 5, ncenters
!!$            do j = 1, 3
!!$               initpts(3,j) = rxyz(j,list(i))
!!$            end do
!!$            call volume_tetrahedron(inipts,volume,stat)
!!$            if(stat == 0) exit  !we have found our point
!!$            if(stat .ne. 0 and i == ncenters) collinear = .true.
!!$         end do        
!!$      end if
!!$      if(ncenters == 4 .or. collinear) then !we have a triangle
!!$      nfacets = 1
!!$      call graham_convex_hull()  !find the extremum points of the line
!!$      do i = 1, 3
!!$         facets(1,i) = list(i)
!!$      end do 
!!$      end if
!!$    else if(stat == 3) then !all four points are collinear
!!$!!      Should never happen so code it another time
!!$!!      if(ncenter >= 7) then
!!$!!      end if
!!$!!      if(ncenter < 6) then
!!$!!        nfacets = 0
!!$!!        return
!!$!!      end if
!!$!!        
!!$!!         if(ncenter == 7) then
!!$!!           
!!$!!      end if
!!$      stop 'All initial points are collinear!'        
!!$   end if
!!$      
!!$
!!$   do ipts =  1, ncenter - 4
!!$   !2) For each remaining points
!!$
!!$      !a) determine the set of facets visible by the point
!!$      !   For this, only need to check the normal of the planes
!!$      
!!$
!!$      !b) determine the set of horizon edges
!!$
!!$      !c) For each horizon edge construct a new triangular facet connecting edge and point p
!!$
!!$      !d) Discard all facets previously visible to p
!!$    end do
!!$
!!$end subroutine convex_hull_construction_3D_CS1989

!!$subroutine volume_tetrahedron(rxyz,volume,stat)
!!$   use module_types
!!$   implicit none
!!$   real(gp), dimension(3,4), intent(in) :: rxyz
!!$   real(gp), intent(out) :: volume
!!$   integer, intent(out) :: stat   ! = 0 volume not zero, = 1 if crossprod is zero, = 2 dot prod gives zero, =3 all points collinear
!!$   !Local variables
!!$   integer :: i
!!$   real(gp), dimension(3) :: vec1, vec2, vec3, crossprod
!!$   real(gp) :: proj1, proj2, proj3
!!$
!!$   ! Initialize stat has all ok
!!$   stat = 0
!!$
!!$   !   volume = (x_3 - x_1) dot [(x_2 - x_1) cross (x_4 - x_3)]
!!$   do i = 1, 3
!!$      vec1(i) = rxyz(i,3) - rxyz(i,1)
!!$      vec2(i) = rxyz(i,2) - rxyz(i,1)
!!$      vec3(i) = rxyz(i,4) - rxyz(i,3)
!!$   end do
!!$
!!$   crossprod(1) = vec2(2)*vec3(3) - vec2(3)*vec3(2)
!!$   crossprod(2) = vec2(3)*vec3(1) - vec2(1)*vec3(3)
!!$   crossprod(3) = vec2(1)*vec3(2) - vec2(2)*vec3(1)
!!$ 
!!$   if(crossprod(1)**2+crossprod(2)**2+crossprod(3)**2 < 1.0d-6) stat = 2
!!$
!!$   volume = 0.0
!!$   do i = 1, 3
!!$      volume =  volume + vec1(i)*crossprod(i)
!!$   end do
!!$
!!$   if(volume < 1.0d-6 .and. stat .ne. 1) stat = 1
!!$
!!$   !Should check if all points are collinear
!!$   proj1 = 0.0
!!$   proj2 = 0.0
!!$   proj3 = 0.0
!!$   do i = 1, 3
!!$      proj1 = proj1 + vec1(i)*vec2(i) / (sqrt(vec1(1)**2+vec1(2)**2+vec1(3)**2)*sqrt(vec2(1)**2+vec2(2)**2+vec2(3)**2))
!!$      proj2 = proj2 + vec1(i)*vec3(i) / (sqrt(vec1(1)**2+vec1(2)**2+vec1(3)**2)*sqrt(vec3(1)**2+vec3(2)**2+vec3(3)**2))
!!$      proj3 = proj3 + vec2(i)*vec3(i) / (sqrt(vec2(1)**2+vec2(2)**2+vec2(3)**2)*sqrt(vec3(1)**2+vec3(2)**2+vec3(3)**2))
!!$   end do
!!$
!!$   if(abs(proj1)-1.0 < 1.0d-6 .and. abs(proj2)-1.0 < 1.0d-6 .and. abs(proj3)-1.0 < 1.0d-6) stat = 3
!!$
!!$end subroutine volume_tetrahedron
!!$
!!$
!!$subroutine build_correct_triangles_for_projection(natom, rxyz, projN, projS, facet,refpos, southref, CM, rad, inverted, pts)
!!$use module_types
!!$implicit none
!!$integer, intent(in) :: natom                              ! Number of atoms 
!!$real(gp),intent(in) :: rad                                ! radius of the sphere for projection
!!$real(gp), dimension(3),intent(in) :: refpos, southref,CM
!!$real(gp), dimension(3, natom), intent(in) :: rxyz         ! Position fo atoms before projection
!!$real(gp), dimension(natom, 3), intent(in) :: projN        ! atom positions in the northern projection
!!$real(gp), dimension(natom, 3), intent(in) :: projS        ! atom positions in the southern projection
!!$integer, dimension(3), intent(in) :: facet                ! atom index of the triangles forming the surface
!!$logical, intent(out) :: inverted
!!$real(gp), dimension(3), intent(out) :: pts
!!$!Local variables
!!$logical :: check
!!$integer :: iface, i, j, l, NeglectPoint
!!$integer, dimension(3,2) :: pairs
!!$real(gp), dimension(3) :: normal
!!$real(gp), dimension(3,3) :: Nnormals, Snormals
!!$real(gp), dimension(4,3) :: New_pts, new_projN, new_projS
!!$
!!$print *,'facet: ',facet
!!$
!!$! Index of the points defining the three edges of the triangles
!!$pairs(1,1) = 1 ; pairs(1,2) = 2
!!$pairs(2,1) = 1 ; pairs(2,2) = 3
!!$pairs(3,1) = 2 ; pairs(3,2) = 3
!!$
!!$! Calculate interior normals for the northern projection
!!$call calculate_interior_normals(natom, projN, facet, Nnormals)
!!$
!!$! Calculate the interior normal of the edges for the southern projection
!!$call calculate_interior_normals(natom, projS, facet, Snormals)  
!!$do i = 1, 3
!!$print *,'ProjN: ',projN(facet(i),:)
!!$print *,'ProjS: ',projS(facet(i),:)
!!$end do
!!$print *,'NORMALS TEST:', dot_product(Nnormals(1,:),Snormals(1,:)) > 0,&
!!$        dot_product(Nnormals(2,:),Snormals(2,:)) > 0 ,&
!!$        dot_product(Nnormals(3,:),Snormals(3,:)) > 0
!!$do i=1,3
!!$print *,'Normals:', Nnormals(i,:),Snormals(i,:)
!!$end do
!!$
!!$! Compare the two, if it does not invert, divide the triangle along this edge 
!!$! They either completly invert, or they do not
!!$! Use greater then zero, because switching poles incures a 180 rotation
!!$if(dot_product(Nnormals(1,:),Snormals(1,:)) > 0 .or. dot_product(Nnormals(2,:),Snormals(2,:)) > 0 .or.&
!!$    dot_product(Nnormals(3,:),Snormals(3,:)) > 0) then
!!$   ! Must select the edge to split
!!$   ! FOR NOW: try splitting edge by edge until we find the one opposite the concave point (the one which builds 2 proprer oriented triangles)
!!$   loop_edge : do i = 1, 3
!!$
!!$         ! Split the triangle with respect to this edge
!!$         do j = 1, 3
!!$            New_pts(1,j) = (rxyz(j,facet(pairs(i,1))) + rxyz(j,facet(pairs(i,2))))/2
!!$            do l = 1, 3
!!$               New_pts(l+1,j) = rxyz(j,facet(l))
!!$            end do
!!$         end do
!!$
!!$         ! Now project the new points
!!$         call stereographic_projection(4, New_pts, refpos, CM, rad, new_projN, normal, NeglectPoint)
!!$         if(NeglectPoint .ne. 0) stop 'Neglecting one of the points'
!!$         call stereographic_projection(4, New_pts, southref, CM, rad, new_projS, normal, NeglectPoint)
!!$         if(NeglectPoint .ne. 0) stop 'Neglecting one of the points'
!!$
!!$         ! Check if the new triangles solved the problem
!!$         ! First triangle
!!$         call calculate_interior_normals(4, new_projN, (/ 1, 2, 4/), Nnormals)
!!$         call calculate_interior_normals(4, new_projS, (/ 1, 2, 4 /), Snormals)
!!$         check = dot_product(Nnormals(i,:),Snormals(i,:)) > 0
!!$         ! Second triangle
!!$         call calculate_interior_normals(4, new_projN, (/ 1, 3, 4/), Nnormals)
!!$         call calculate_interior_normals(4, new_projS, (/ 1, 3, 4 /), Snormals)
!!$         check = check .and. dot_product(Nnormals(i,:),Snormals(i,:)) > 0
!!$
!!$         ! If the triangles are now well behaving, keep this configuration
!!$         if(check) then
!!$            inverted = .true.
!!$            do j = 1, 3 
!!$               pts(j) = new_projN(1,j)
!!$            end do
!!$            exit loop_edge
!!$         end if
!!$         if(.not. check .and. i==3) stop 'Could not find a correct convex representation for the stereographic graph!'
!!$   end do loop_edge
!!$
!!$else
!!$   inverted = .false.
!!$   pts = (/ 1.0, 1.0, 1.0 /)      
!!$end if
!!$
!!$end subroutine build_correct_triangles_for_projection
!!$
!!$subroutine calculate_interior_normals(natom, proj, facet, normals)
!!$use module_types
!!$implicit none
!!$integer, intent(in) :: natom                         ! Number of atoms 
!!$real(gp), dimension(natom, 3), intent(in) :: proj    ! atom positions in the projection
!!$integer, dimension(3), intent(in) :: facet           ! atom index of the triangle forming the surface
!!$real(gp), dimension(3,3), intent(out) :: normals
!!$! Local variables
!!$integer :: i, j, k, l, pt1, pt2
!!$real(gp), dimension(3) :: norm
!!$real(gp), dimension(3,3) :: edges
!!$print *,'cin:facet',facet
!!$   ! Calculate the edges for the northern projection
!!$   l = 0
!!$   do i = 1, 3
!!$      do j = i, 3
!!$         if (i == j) cycle
!!$         l = l + 1
!!$         pt1 = facet(i)
!!$         pt2 = facet(j)
!!$         do k = 1, 3
!!$            edges(l,k) = proj(pt2,k)-proj(pt1,k)   !numbering of the edges 1=1-2, 2=1-3, 3=2-3 
!!$         end do
!!$      end do
!!$   end do
!!$
!!$   !calculate norms
!!$   do i =1, 3
!!$     do j = 1,3
!!$        norm(i) = norm(i) + edges(i,j)**2 
!!$     end do
!!$   end do
!!$
!!$   ! The interior normals to these edges are defined by
!!$   do j = 1, 3
!!$      normals(1,j) = proj(facet(3),j) - edges(1,j)*dot_product(edges(2,:),edges(1,:))&!(edges(2,1)*edges(1,1)+edges(2,2)*edges(1,2)+edges(2,3)*edges(1,3))&
!!$                     / norm(1) - proj(facet(1),j)
!!$      normals(2,j) = proj(facet(2),j) - edges(2,j)*(edges(2,1)*edges(1,1)+edges(2,2)*edges(1,2)+edges(2,3)*edges(1,3))&
!!$                     / norm(2) - proj(facet(1),j)
!!$      normals(3,j) = proj(facet(1),j) - edges(3,j)*(-edges(1,1)*edges(3,1)-edges(1,2)*edges(3,2)-edges(1,3)*edges(3,3))&
!!$                     / norm(3) - proj(facet(2),j) 
!!$   end do
!!$
!!$end subroutine calculate_interior_normals

