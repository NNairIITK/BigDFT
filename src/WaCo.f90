!> @file
!! Wannier constructor
!! @author
!!    Copyright (C) 2011-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Program to construct the Wannier functions
program WaCo

   use module_base
   use module_types
   use module_interfaces, except_this_one => writeonewave
   use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
   use yaml_output
   use module_input_dicts
   use module_atoms, only: deallocate_atoms_data
   use communications_base, only: comms_cubic
   use communications_init, only: orbitals_communicators
   implicit none
   character :: filetype*4,outputype*4
   type(locreg_descriptors) :: Glr
   type(orbitals_data) :: orbs,orbsw,orbsv
   type(atoms_data) :: atoms
   type(input_variables) :: input
   type(workarr_sumrho) :: w
   type(comms_cubic), target :: commsw
   type(local_zone_descriptors) :: Lzd             !< debug only
   integer :: iiband,ldim,gdim                     !< debug only
   logical, dimension(:),allocatable :: calcbounds !< debug only
   real(gp), parameter :: b2a=0.5291772108_dp
   real :: tcpu0,tcpu1
   real(gp) :: tel
   real(gp), dimension(3) :: shift,CM
   real(gp) :: dist,rad,sprdfact,sprddiff,enediff,sprdmult
   integer :: iproc, nproc, i_stat, ierr, npsidim
   integer :: nvirtu,nvirtd,nrpts
   integer :: NeglectPoint, CNeglectPoint
   integer :: ncount0,ncount1,ncount_rate,ncount_max,iat,iformat
   integer :: iwann, iiwann, iband, nwann, nband, plotwann
   integer :: iw1, iw2, ifile
   integer :: nsprd,ndiag, nwannCon
   character(len=*), parameter :: subname='WaCo'
   character(len=4) :: num, units
   integer, allocatable :: wann_list(:), Zatoms(:,:), ncenters(:), wtypes(:,:)
   real(gp), dimension(:,:), pointer :: rxyz_old, cxyz,rxyz_wann
   real(gp), dimension(:,:), allocatable :: radii_cf
   real(gp), allocatable :: sprd(:), locrad(:), eigen(:,:), proj(:,:), projC(:,:),distw(:),charge(:),prodw(:),wannocc(:)
   real(wp), allocatable :: psi(:,:),wann(:),wannr(:),lwann(:)
   real(wp), allocatable :: ham(:,:,:),hamr(:,:,:)
   real(wp), allocatable :: diag(:,:),diagT(:)
   integer, dimension(:), pointer :: buf
   character(len=60) :: radical, filename, run_id
   logical :: notocc, bondAna,Stereo,hamilAna,WannCon,linear,outformat
   integer, dimension(:), allocatable :: ConstList
   integer, allocatable :: nfacets(:),facets(:,:,:),vertex(:,:,:), l(:), mr(:)
   real(gp), dimension(3) :: refpos, normal, box
   real(kind=8),dimension(:,:),allocatable :: umn, umnt, rho, rhoprime, amn, tmatrix
   integer :: i, j, k, i_all, ilr
   character(len=16) :: seedname
   integer :: n_occ, n_virt, n_virt_tot, nproj,nband_old,nkpt_old,iwann_out 
   logical :: w_unk, w_sph, w_ang, w_rad, pre_check,residentity,write_resid
   integer, allocatable, dimension (:) :: virt_list
   integer :: nbandCon,nconfig
   integer, dimension(:),allocatable :: bandlist
   logical :: idemp
   integer, dimension(4) :: mpi_info
   type(dictionary), pointer :: user_inputs
   external :: gather_timings
   ! ONLY FOR DEBUG
!   real(gp) :: Gnorm, Lnorm
!   integer :: indL,ilr
!   real(kind=8),dimension(:,:),allocatable :: coeff
!   type(orbitals_data) :: wannorbs
!   type(comms_cubic), target :: wanncomms
!   real(wp), allocatable :: psi2(:)
!   real(gp),allocatable :: cxyz2(:,:) !debug only
!   integer, allocatable :: list(:)    !debug only
   interface
      subroutine read_inter_header(iproc,seedname, filetype,residentity,write_resid, n_occ, pre_check,&
                 n_virt_tot, n_virt, w_unk, w_sph, w_ang, w_rad)
        implicit none
        integer, intent(in) :: iproc
        character, intent(out) :: seedname*16, filetype*4
        integer, intent(out) :: n_occ, n_virt, n_virt_tot
        logical, intent(out) :: w_unk, w_sph, w_ang, w_rad, pre_check,residentity,write_resid
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

   !-finds the number of taskgroup size
   !-initializes the mpi_environment for each group
   !-decides the radical name for each run
   call bigdft_init(mpi_info,nconfig,run_id,ierr)

   !just for backward compatibility
   iproc=mpi_info(1)
   nproc=mpi_info(2)

   if (nconfig < 0) stop 'runs-file not supported for WaCo executable'

   call dict_init(user_inputs)
   call user_dict_from_files(user_inputs, trim(run_id)//trim(bigdft_run_id_toa()), &
        & 'posinp'//trim(bigdft_run_id_toa()), bigdft_mpi)
   call inputs_from_dict(input, atoms, user_inputs)
   call dict_free(user_inputs)

!!$   if (input%verbosity > 2) then
!!$      nproctiming=-nproc !timing in debug mode
!!$   else
!!$      nproctiming=nproc
!!$   end if

   !call timing(nproctiming,'WaCo_time.prc','IN')
   call f_timing_reset(filename=trim(input%dir_output)//'WaCo_time.yaml',&
        master=iproc==0,&
        verbose_mode=input%verbosity>2)


   call cpu_time(tcpu0)
   call system_clock(ncount0,ncount_rate,ncount_max) 
   call timing(iproc,'Precondition  ','ON')   

   !###################################################################
   ! Read input files and initialise the variables for the wavefunctions
   !###################################################################
!   call standard_inputfile_names(input,radical,nproc)

   call Waco_input_variables(iproc,trim(radical)//'.waco',nband,nwann,bondAna,Stereo,hamilAna,WannCon,&
        outputype,nwannCon,refpos,units,sprdfact,sprddiff,enediff,outformat,linear,nbandCon,sprdmult)

   if(.not.linear) nbandCon=1
   allocate(ConstList(nwannCon),stat=i_stat)
   call memocc(i_stat, ConstList,'ConstList',subname)
   allocate(bandlist(nbandCon),stat=i_stat)
   call memocc(i_stat, bandlist,'bandlist',subname)

   call read_input_waco(trim(radical)//'.waco',nwannCon,ConstList,linear,nbandCon,bandlist) 

   allocate(radii_cf(atoms%astruct%ntypes,3+ndebug),stat=i_stat)
   call memocc(i_stat,radii_cf,'radii_cf',subname)

   call system_properties(iproc,nproc,input,atoms,orbs,radii_cf)

   ! Determine size alat of overall simulation cell and shift atom positions
   ! then calculate the size in units of the grid space
   call system_size(atoms,atoms%astruct%rxyz,radii_cf,input%crmult,input%frmult,input%hx,input%hy,input%hz,&
        .false.,Glr,shift)
   if (iproc == 0) &
        & call print_atoms_and_grid(Glr, atoms, atoms%astruct%rxyz, shift, input%hx,input%hy,input%hz)
   
   box(1) = atoms%astruct%cell_dim(1)*b2a !Glr%d%n1*input%hx * b2a
   box(2) = atoms%astruct%cell_dim(2)*b2a !Glr%d%n2*input%hy * b2a
   box(3) = atoms%astruct%cell_dim(3)*b2a !Glr%d%n3*input%hz * b2a

   ! Create wavefunctions descriptors and allocate them inside the global locreg desc.
   call createWavefunctionsDescriptors(iproc,input%hx,input%hy,input%hz,&
        atoms,atoms%astruct%rxyz,radii_cf,input%crmult,input%frmult,Glr)
   if (iproc == 0) call print_wfd(Glr%wfd)

   ! don't need radii_cf anymore
   i_all = -product(shape(radii_cf))*kind(radii_cf)
   deallocate(radii_cf,stat=i_stat)
   call memocc(i_stat,i_all,'radii_cf',subname)

   !#################################################################
   ! Read Other files
   !#################################################################
   !input.inter
   call read_inter_header(iproc,seedname, filetype,residentity,write_resid, n_occ, pre_check,&
        n_virt_tot, n_virt, w_unk, w_sph, w_ang, w_rad)
   allocate(virt_list(n_virt),stat=i_stat)
   call memocc(i_stat,virt_list,'virt_list',subname)
   if (n_virt .ne. 0) then
      call read_inter_list(iproc, n_virt, virt_list)
   end if 

   if ((nband .ne. n_occ+n_virt) .and. iproc == 0) then
      call yaml_warning('Number of bands in the .waco file' // trim(yaml_toa(nband)))
      call yaml_warning('not equal to the number of bands used in .inter file' // trim(yaml_toa(n_occ+n_virt)))
      !write(*,*) 'Number of bands in the .waco file : ',nband
      !write(*,*) 'not equal to the number of bands used in .inter file:', n_occ+n_virt
      call mpi_finalize(ierr)
      stop
   end if

   ! seedname.umn
   allocate(umn(nwann,nband),stat=i_stat)
   call memocc(i_stat,umn,'umn',subname)
   call read_umn(iproc,nwann,nband,seedname,umn)

   !seedname.dat
   allocate(wann_list(nwann),stat=i_stat)
   call memocc(i_stat,wann_list,'wann_list',subname)
   call read_spread_file(iproc,seedname,nwann,plotwann,wann_list)

   ! Read Wannier centers
   allocate(cxyz(3,plotwann),stat=i_stat)
   call memocc(i_stat,cxyz,'cxyz',subname)
   allocate(rxyz_wann(3,atoms%astruct%nat),stat=i_stat)
   call memocc(i_stat,rxyz_wann,'rxyz_wann',subname)
   allocate(sprd(plotwann+1),stat=i_stat)
   call memocc(i_stat,sprd,'sprd',subname)
   call read_centers(iproc,nwann,plotwann,atoms%astruct%nat,seedname,wann_list,cxyz,rxyz_wann,.true.,sprd)

   call timing(iproc,'Precondition  ','OF')

   !###################################################################
   ! Constructing the density matrix
   !###################################################################

   if (iproc == 0) then
      call yaml_comment('Constructing density matrix',hfill='-')
      !write(*,*) '!==================================!'
      !write(*,*) '!   Constructing density matrix    !'
      !write(*,*) '!==================================!'
      !write(*,*)
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

!   if (iproc == 0) then
!      do i=1,nwann
!         do j=1,nwann
!            write(*,*) i , j ,rho(i,j)
!         end do
!      end do
!   end if

   if (iproc == 0) then
      call yaml_map('Constructing density matrix',.true.)
      !write(*,*) '!===================================!'
      !write(*,*) '!Constructing density matrix : DONE !'
      !write(*,*) '!===================================!'
      !write(*,*)
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
           call yaml_warning('Not indempotent (' // trim(yaml_toa((/i,j/))) &
                & // trim(yaml_toa( (/ rhoprime(i,j), rho(i,j) /))) // ')')
           !write(*,*) 'Not indempotent',i,j,rhoprime(i,j), rho(i,j)
           idemp = .false.
         end if
      end do
   end do

   if(idemp .eqv. .false.) then
      call yaml_warning('Density matrix not idempotent')
     stop 
   end if

   i_all = -product(shape(rhoprime))*kind(rhoprime)
   deallocate(rhoprime,stat=i_stat)
   call memocc(i_stat,i_all,'rhoprime',subname)
   i_all = -product(shape(rho))*kind(rho)
   deallocate(rho,stat=i_stat)
   call memocc(i_stat,i_all,'rho',subname)


   !########################################################
   ! Bonding analysis
   !########################################################
  if(bondAna) then
     if (iproc == 0) then
         call yaml_comment('Bonding Analysis',hfill='-')
         call yaml_open_map('Bonding Analysis')
         !write(*,*) '!==================================!'
         !write(*,*) '!     Bonding Analysis :           !'
         !write(*,*) '!==================================!'
     end if


      !wann_list will contain the list of occupied Wannier functions
      allocate(Zatoms(atoms%astruct%nat,plotwann),stat=i_stat)
      call memocc(i_stat,Zatoms,'Zatoms',subname)  
      allocate(ncenters(plotwann),stat=i_stat)
      call memocc(i_stat,ncenters,'ncenters',subname)


      ! Now calculate the bonding distances and ncenters
      if (iproc == 0) then
         call yaml_open_sequence('Number of atoms associated to the WFs')
         call yaml_comment('WF   OWF   Nc   Spr(ang^2)   Atom numbers')
         !write(*,'(2x,A)') 'Number of atoms associated to the WFs'
         !write(*,'(3x,A,1x,A,4x,A,3x,A,3x,A)') 'WF','OWF','Spr(ang^2)','Nc','Atom numbers:'
      end if
      Zatoms = 0
      do iwann = 1, plotwann
         ncenters(iwann) = 0
         iat = 0
         do i = 1, atoms%astruct%nat
            call get_mindist(Glr%geocode,rxyz_wann(1,i),cxyz(1,iwann),box,dist)
            if (dist**2 <= sprdfact * sprd(iwann)) then    !for normal distribution: 1=68%, 1.64=80%, 3=94%
               ncenters(iwann) = ncenters(iwann) +1
               iat = iat +1
               Zatoms(iat,iwann) = i
            end if
         end do
         if(iproc == 0) then
           call yaml_sequence(advance="no")
           call yaml_open_sequence(flow=.true.)
              call yaml_sequence(trim(yaml_toa( (/ wann_list(iwann), iwann, ncenters(iwann) /) )))
              call yaml_sequence(trim(yaml_toa(sprd(iwann),fmt='(f14.6)')))
              call yaml_sequence(trim(yaml_toa(Zatoms(1:iat,iwann))))
           call yaml_close_sequence()
           !write(*,'(2I4,F14.6,2x,I4,6(2x,I4))') &
           !   & wann_list(iwann),iwann, sprd(iwann), ncenters(iwann), (Zatoms(i_all,iwann),i_all=1,iat)
         end if
      end do
      if (iproc == 0) call yaml_close_sequence()

     ! Calculate occupation of the wannier functions
     allocate(wannocc(nwann),stat=i_stat)
     call memocc(i_stat,wannocc,'wannocc',subname)
     wannocc = 0.0_dp
     do iwann = 1, nwann
        do iband = 1, n_occ
            wannocc(iwann) = wannocc(iwann) + umn(iwann,iband)**2
        end do
     end do
     if(iproc == 0) call yaml_map('Total number of electrons',2*sum(wannocc))
     !if(iproc == 0) write(*,*) 'Total number of electrons: ',2*sum(wannocc)

      allocate(distw(maxval(ncenters)),stat=i_stat)
      call memocc(i_stat,distw,'distw',subname)
      allocate(charge(atoms%astruct%nat),stat=i_stat)
      call memocc(i_stat,charge,'charge',subname)
      allocate(prodw(maxval(ncenters)),stat=i_stat)
      call memocc(i_stat,prodw,'prodw',subname)

      ! Calculate "Wannier charge" = 2 *occupation* product(d_all_atoms__associated_with_this_wannier_except_this_one) / d_total
      charge = 0.0_dp
      do iwann = 1, plotwann
         distw = 0.0_dp
         do iat = 1, ncenters(iwann)
            call get_mindist(Glr%geocode,rxyz_wann(1,Zatoms(iat,iwann)),cxyz(1,iwann),box,distw(iat))
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

!      do iat = 1, atoms%astruct%nat
!         write(*,*) 'Charge of atom ', iat,' is : ', charge(iat)
!      end do
      if(iproc == 0) then
         call yaml_map('Total Wannier charge',sum(charge))
         call yaml_open_map('Maximal charge',flow=.true.)
            call yaml_map('Atom', maxloc(charge))
            call yaml_map('Value',maxval(charge))
         call yaml_close_map()
         call yaml_open_map('Minimal charge',flow=.true.)
            call yaml_map('Atom', minloc(charge))
            call yaml_map('Value',minval(charge))
         call yaml_close_map()
         !write(*,*) 'Total Wannier charge is : ',sum(charge)
         !write(*,*) 'Maximal charge is: ',maxval(charge), 'on atom :', maxloc(charge)
         !write(*,*) 'Minimal charge is: ',minval(charge), 'on atom :', minloc(charge)  
      end if

      open(unit=22,file='Wannier_charge.dat',status='unknown')
      do iwann = 1, plotwann
          write(22, '(E14.6, 2x, E14.6)') 0.0_dp, 0.0_dp
      end do
      do iat = 1, atoms%astruct%nat
         write(22, '(E14.6, 2x, E14.6)') (charge(iat)-minval(charge))/(maxval(charge)-minval(charge)), charge(iat)
      end do
      close(unit=22)

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

      if(iproc == 0) then
         call character_list(nwann,nproj,tmatrix,plotwann,ncenters,wann_list,l,mr) 
      end if
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
!     write(*,*) '######################ENTERING DEBUG#################'
!      allocate(distw(atoms%astruct%nat),stat=i_stat)
!      call memocc(i_stat,distw,'distw',subname)
!      allocate(charge(atoms%astruct%nat),stat=i_stat)
!      call memocc(i_stat,charge,'charge',subname)
!      allocate(prodw(atoms%astruct%nat),stat=i_stat)
!      call memocc(i_stat,prodw,'prodw',subname)
!      allocate(cxyz2(3,nwann),stat=i_stat)
!      call memocc(i_stat,cxyz2,'cxyz2',subname)
!      allocate(list(nwann),stat=i_stat)
!      call memocc(i_stat,list,'list',subname)
!      do iwann = 1, nwann
!         list(iwann) = iwann
!      end do 
!      call read_centers(iproc,nwann,nwann,atoms%astruct%nat,seedname,list,cxyz2,rxyz_wann,.false.,sprd)
!      charge = 0.0_dp
!      do iwann = 1, nwann
!         distw = 0.0_dp
!         do iat = 1, atoms%astruct%nat
!            distw(iat) = (rxyz_wann(1,iat)-cxyz2(1,iwann))**2 +&
!                        (rxyz_wann(2,iat)-cxyz2(2,iwann))**2 +&
!                        (rxyz_wann(3,iat)-cxyz2(3,iwann))**2
!         end do
!         prodw = 0.0_dp
!         do iat = 1, atoms%astruct%nat
!            prodw(iat) = 1.0_dp
!            do i = 1, atoms%astruct%nat
!               if(i == iat) cycle
!write(*,*) 'distw(i),prodw(iat)',distw(i),prodw(iat)
!               prodw(iat) = prodw(iat) * distw(i)
!            end do
!         end do
!         do iat = 1, atoms%astruct%nat
!            charge(iat) = charge(iat) + 2.0_dp * wannocc(iwann) * prodw(iat)  / (sum(prodw))
!         end do
!      end do
!      if (iproc == 0) then
!         write(*,*) 'Total Wannier(all) charge is : ',sum(charge)
!         write(*,*) 'Maximal charge is: ',maxval(charge), 'on atom :', maxloc(charge)
!         write(*,*) 'Minimal charge is: ',minval(charge), 'on atom :', minloc(charge)
!      end if
!      open(unit=22,file='Wannier_charge2.dat',status='unknown')
!      do iat = 1, atoms%astruct%nat
!         write(22, '(E14.6, 2x, E14.6)') charge(iat)/maxval(charge), charge(iat)
!      end do
!      close(unit=22)
!      i_all = -product(shape(distw))*kind(distw)
!      deallocate(distw,stat=i_stat)
!      call memocc(i_stat,i_all,'distw',subname)   
!      i_all = -product(shape(charge))*kind(charge)
!      deallocate(charge,stat=i_stat)
!      call memocc(i_stat,i_all,'charge',subname)   
!      i_all = -product(shape(prodw))*kind(prodw)
!      deallocate(prodw,stat=i_stat)
!      call memocc(i_stat,i_all,'prodw',subname)   
!END DEBUG

      i_all = -product(shape(wannocc))*kind(wannocc)
      deallocate(wannocc,stat=i_stat)
      call memocc(i_stat,i_all,'wannocc',subname)   

      !DOS of the sprd
      allocate(wtypes(2,nwann))
      wtypes = 1
      call wannier_dos('sprd.dat',1,plotwann,2,wtypes,sprd)
      deallocate(wtypes)

!      call scalar_kmeans_diffIG(0,maxval(sprd(1:plotwann-1))*1.0d-1,plotwann,sprd,'spread',nsprd,buf)
      call scalar_kmeans_diffIG(iproc,0,sprddiff,plotwann,sprd,'spread',nsprd,buf)

!###############################################################
! Stereographic projection stuff
!###############################################################
      if(Stereo) then
         ! Transform rxyz_wann from angstrom to bohr
         do iat = 1, atoms%astruct%nat
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

         ! Output the posref for visual check
         open(unit=22,file='pos_ref.xyz', status='unknown')
         write(22,'(I4)') atoms%astruct%nat+1
         write(22,*) !skip this line
         write(22,'(A,3(2x,E14.6))') 'X',(refpos(i),i=1,3)
         do i = 1, atoms%astruct%nat
            write(22,'(A,3(2x,E14.6))')atoms%astruct%atomnames(atoms%astruct%iatype(i)),rxyz_wann(1,i),rxyz_wann(2,i),rxyz_wann(3,i)
         end do
         close(unit=22)

         ! Calculate Center of mass
         CM = 0.0_dp
         do iat = 1, atoms%astruct%nat
            do j = 1, 3 
               CM(j) = CM(j) + rxyz_wann(j,iat) / real(atoms%astruct%nat,kind=8)
            end do
         end do

         !Calculate the radius of the sphere (choose it to be the biggest distance from the CM)
         rad= 0.0_dp
         do iat = 1, atoms%astruct%nat
            dist = 0.0_dp
            do j = 1, 3
               dist = dist + (rxyz_wann(j,iat) - CM(j))**2
            end do
            rad = max(dist, rad)
         end do
         rad =sqrt(rad)

         allocate(proj(atoms%astruct%nat,3),stat=i_stat)
         call memocc(i_stat,proj,'proj',subname)
         allocate(projC(plotwann,3),stat=i_stat)
         call memocc(i_stat,projC,'projC',subname)
         allocate(nfacets(plotwann),stat=i_stat)
         call memocc(i_stat,nfacets,'nfacets',subname)
         allocate(facets(plotwann,maxval(ncenters)*(maxval(ncenters)-1)/2,3),stat=i_stat)
         call memocc(i_stat,facets,'facets',subname)
         allocate(vertex(plotwann,maxval(ncenters)*(maxval(ncenters)-1)/2,3),stat=i_stat)
         call memocc(i_stat,vertex,'vertex',subname)

         ! Do stereographic projection of atoms and Wannier centers
         call stereographic_projection(0,atoms%astruct%nat,rxyz_wann,refpos, CM, rad, proj, normal, NeglectPoint)
         ! TO DO: CNeglectPoint should be a vector...
         call stereographic_projection(1,plotwann,cxyz,refpos, CM, rad, projC, normal, CNeglectPoint)
         call shift_stereographic_projection(plotwann,projC,atoms%astruct%nat,proj)
         call write_stereographic_projection(22, 'proj.xyz    ', atoms, proj, NeglectPoint) 

         !Must warn if a Wannier center is on the projection reference
         if(CNeglectPoint .ne. 0) then
            call yaml_warning('The Wannier center' // trim(yaml_toa(CNeglectPoint)) // &
               & ' is on the refence point of the projection.')
            call yaml_warning('Surfaces will be deformed.')
            !write(*,*) 'The Wannier center ',CNeglectPoint,'is on the refence point of the projection.'
            !write(*,*) 'Surfaces will be deformed'
!            call mpi_finalize(ierr)
!            stop
         end if
         
         call build_stereographic_graph_facets(atoms%astruct%nat,plotwann,maxval(ncenters),4.0d0,rxyz_wann,ncenters,&
              Zatoms,nfacets,facets,vertex)
         call output_stereographic_graph(atoms%astruct%nat,maxval(ncenters),proj,projC,plotwann,ncenters,Zatoms,nfacets,&
              facets,vertex,normal,NeglectPoint)

         i_all = -product(shape(nfacets))*kind(nfacets)
         deallocate(nfacets,stat=i_stat)
         call memocc(i_stat,i_all,'nfacets',subname)
         i_all = -product(shape(facets))*kind(facets)
         deallocate(facets,stat=i_stat)
         call memocc(i_stat,i_all,'facets',subname)
         i_all = -product(shape(vertex))*kind(vertex)
         deallocate(vertex,stat=i_stat)
         call memocc(i_stat,i_all,'vertex',subname)
         i_all = -product(shape(proj))*kind(proj)
         deallocate(proj,stat=i_stat)
         call memocc(i_stat,i_all,'proj',subname)
         i_all = -product(shape(projC))*kind(projC)
         deallocate(projC,stat=i_stat)
         call memocc(i_stat,i_all,'projC',subname)
      end if
      i_all = -product(shape(cxyz))*kind(cxyz)
      deallocate(cxyz,stat=i_stat)
      call memocc(i_stat,i_all,'cxyz',subname)
      i_all = -product(shape(rxyz_wann))*kind(rxyz_wann)
      deallocate(rxyz_wann,stat=i_stat)
      call memocc(i_stat,i_all,'rxyz_wann',subname)
      i_all = -product(shape(sprd))*kind(sprd)
      deallocate(sprd,stat=i_stat)
      call memocc(i_stat,i_all,'sprd',subname)
      i_all = -product(shape(Zatoms))*kind(Zatoms)
      deallocate(Zatoms,stat=i_stat)
      call memocc(i_stat,i_all,'Zatoms',subname)
      if(.not. hamilAna) then
         i_all = -product(shape(ncenters))*kind(ncenters)
         deallocate(ncenters,stat=i_stat)
         call memocc(i_stat,i_all,'ncenters',subname)
      end if

     if (iproc == 0) then
         call yaml_close_map()
         !write(*,*) '!==================================!'
         !write(*,*) '!     Bonding Analysis : DONE      !' 
         !write(*,*) '!==================================!'
     end if

  end if
  if (hamilAna) then
     if (iproc == 0) then
        call yaml_comment('Hamiltonian Analysis',hfill='-') 
        call yaml_open_map('Hamiltonian Analysis') 
        !write(*,*) '!==================================!'
        !write(*,*) '!     Hamiltonian Analysis :       !' 
        !write(*,*) '!==================================!'
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
     allocate(wtypes(nrpts,nwann),stat=i_stat)
     call memocc(i_stat,wtypes,'wtypes',subname)
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
           wtypes(i,iwann) = iw1
        end do
     end do

     allocate(eigen(1,nband),stat=i_stat)
     call memocc(i_stat,eigen,'eigen',subname)
     call read_eigenvalues(trim(seedname)//'.eig',nband,1,eigen)
     call wannier_projected_dos('Wannier_projected_dos.dat',nrpts,nwann,nband,umn,nsprd+1,wtypes,eigen)
     call wannier_dos('Wannier_dos.dat',nrpts,nwann,nsprd+1,wtypes,diag)
     i_all = -product(shape(eigen))*kind(eigen)
     deallocate(eigen,stat=i_stat)
     call memocc(i_stat,i_all,'eigen',subname)

     allocate(diagT(plotwann),stat=i_stat)
     call memocc(i_stat,diagT,'diagT',subname)

     if (iproc == 0) call yaml_open_sequence('Diagonal of the Hamiltonian (without empty WFs)')
     !if(iproc == 0) write(*,'(A)') 'Diagonal of the Hamiltonian (without empty WFs)'
     do i = 1, nrpts
        if (iproc == 0) then
           call yaml_sequence(advance='no')
           call yaml_comment(trim(yaml_toa(i,fmt='(i4.4)')))
           call yaml_open_sequence(flow=.true.)
        end if
        do iwann = 1, plotwann
           if (iproc == 0) then
              call yaml_open_sequence(flow=.true.)
              call yaml_sequence(trim(yaml_toa(ham(i,iwann,iwann),fmt='(e14.6')))
              if (bondAna) call yaml_sequence(trim(yaml_toa(ncenters(iwann))))
              call yaml_close_sequence()
              call yaml_comment(trim(yaml_toa(iwann,fmt='(i4.4)')))
           end if
           !if (iproc == 0 .and. bondAna) then
           !   write(*,'(i4,2x,i4,2x,E14.6,2x,i4)') iwann,iwann,ham(i,iwann,iwann),ncenters(iwann)
           !else if (iproc == 0) then
           !   write(*,'(i4,2x,i4,2x,E14.6)') iwann,iwann,ham(i,iwann,iwann)
           !end if
           diagT(iwann) = hamr(i,iwann,iwann)
        end do
        if (iproc == 0) call yaml_close_sequence()
     end do
     if (iproc == 0) call yaml_close_sequence()

     !Deallocate buf, because it is allocated again in scalar_kmeans_diffIG
     i_all=-product(shape(buf))*kind(buf)
     deallocate(buf,stat=i_stat)
     call memocc(i_stat,i_all,'buf',subname)
     if(.not. bondAna) nsprd = 0  !if we didn't do bonding analysis, find here the best for the hamiltonian
     call scalar_kmeans_diffIG(iproc,nsprd,enediff,plotwann,diagT,'diagonal',ndiag,buf)

     i_all=-product(shape(wtypes))*kind(wtypes)
     deallocate(wtypes,stat=i_stat)
     call memocc(i_stat,i_all,'wtypes',subname)
     i_all=-product(shape(diagT))*kind(diagT)
     deallocate(diagT,stat=i_stat)
     call memocc(i_stat,i_all,'diagT',subname)
     i_all=-product(shape(hamr))*kind(hamr)
     deallocate(hamr,stat=i_stat)
     call memocc(i_stat,i_all,'hamr',subname)
     i_all=-product(shape(diag))*kind(diag)
     deallocate(diag,stat=i_stat)
     call memocc(i_stat,i_all,'diag',subname)
     if(bondAna)then
        i_all = -product(shape(ncenters))*kind(ncenters)
        deallocate(ncenters,stat=i_stat)
        call memocc(i_stat,i_all,'ncenters',subname)
     end if

     if (iproc == 0) then 
        call yaml_close_map()
        !write(*,*) '!==================================!'
        !write(*,*) '!     Hamiltonian Analysis : DONE  !' 
        !write(*,*) '!==================================!'
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
         & orbs%nspin,orbs%nspinor,orbs%nkpts,orbs%kpts,orbs%kwgts,orbsw,.false.)
     call orbitals_communicators(iproc,nproc,Glr,orbsw,commsw)

     nvirtu = n_virt
     nvirtd = 0
     if (input%nspin==2) nvirtd=0!nvirtu
     call orbitals_descriptors(iproc,nproc,nvirtu+nvirtd,nvirtu,nvirtd, &
         & orbs%nspin,orbs%nspinor,orbs%nkpts,orbs%kpts,orbs%kwgts,orbsv,.false.)

     if(linear)then
       allocate(cxyz(3,nwannCon),stat=i_stat)
       call memocc(i_stat,cxyz,'cxyz',subname)
       allocate(rxyz_wann(3,atoms%astruct%nat),stat=i_stat)
       call memocc(i_stat,rxyz_wann,'rxyz_wann',subname)
       allocate(sprd(nwannCon+1),stat=i_stat)
       call memocc(i_stat,sprd,'sprd',subname)
       allocate(locrad(nwannCon),stat=i_stat)
       call memocc(i_stat,locrad,'locrad',subname)

       call read_centers(iproc,nwann,nwannCon,atoms%astruct%nat,seedname,ConstList,cxyz,rxyz_wann,.true.,sprd)

       call yaml_open_sequence('Wannier centers')
       do iwann = 1, nwannCon
          locrad(iwann) = sprdmult*(sprd(iwann)/(b2a*b2a))
          if (iproc == 0) then
             call yaml_sequence(advance='no')
             call yaml_open_map('Wannier',flow=.true.)
                call yaml_map('Center',cxyz(1:3,iwann))
                call yaml_map('Locrad',locrad(iwann))
             call yaml_close_map()
             call yaml_comment(trim(yaml_toa(iwann,fmt='(i4.4)')))
             !write(*,*) 'Wannier center for iwann ',iwann,'  :  ',(cxyz(i,iwann),i=1,3)
             !write(*,*) 'Locrad of iwann ',iwann, '  :  ', locrad(iwann)
          end if
       end do
       if (iproc == 0) call yaml_close_sequence

       i_all = -product(shape(rxyz_wann))*kind(rxyz_wann)
       deallocate(rxyz_wann,stat=i_stat)
       call memocc(i_stat,i_all,'rxyz_wann',subname)
       i_all = -product(shape(sprd))*kind(sprd)
       deallocate(sprd,stat=i_stat)
       call memocc(i_stat,i_all,'sprd',subname)
 
     ! Should construct a proper Lzd for each Wannier, then use global -> local transformation
       call nullify_local_zone_descriptors(Lzd)
       call copy_locreg_descriptors(Glr, Lzd%Glr)
       lzd%hgrids(1)=input%hx
       lzd%hgrids(2)=input%hy
       lzd%hgrids(3)=input%hz
       allocate(Lzd%Llr(nwannCon))
       allocate(calcbounds(nwannCon),stat=i_stat)
       call memocc(i_stat,calcbounds,'calcbounds',subname)
       calcbounds =.false.  
       do ilr=1,nwannCon
          Lzd%llr(ilr)%locregCenter(1)=cxyz(1,ilr)
          Lzd%llr(ilr)%locregCenter(2)=cxyz(2,ilr)
          Lzd%llr(ilr)%locregCenter(3)=cxyz(3,ilr)

          Lzd%llr(ilr)%locrad=locrad(ilr)
       end do
       call determine_locregSphere_parallel(iproc,nproc,nwannCon,Lzd%hgrids(1),&
               Lzd%hgrids(2),Lzd%hgrids(3),atoms%astruct,orbs,Lzd%Glr,Lzd%Llr,calcbounds) 
     end if


     !##########################################################################
     ! Read the Wavefunctions
     !##########################################################################

     ! Wavefunctions calculated by BigDFT already are normalized.
     !if (iproc==0) then
     !   write(*,*) '!==================================!'
     !   write(*,*) '!  Reading the wavefunctions       !'
     !   write(*,*) '!==================================!'
     !end if

     call timing(iproc,'CrtProjectors ','ON')

     ! assign the input_wf_format
     iformat = WF_FORMAT_NONE
     select case (filetype)
     case ("ETSF","etsf")
        iformat = WF_FORMAT_ETSF
     case ("BIN","bin")
        iformat = WF_FORMAT_BINARY
     case default
        !if (iproc == 0) write(*,*)' WARNING: Missing specification of wavefunction files'
        call yaml_warning('Missing specification of wavefunction files')
        stop
     end select

     npsidim=max((Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f)*orbsw%norbp*orbsw%nspinor,sum(commsw%ncntt(0:nproc-1)))
     allocate(psi(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,orbsw%norbp*orbsw%nspinor),stat=i_stat)
     call memocc(i_stat,psi,'psi',subname)
     allocate(rxyz_old(3,atoms%astruct%nat),stat=i_stat)  
     call memocc(i_stat,rxyz_old,'rxyz_old',subname)
     ! For the occupied orbitals, need to modifify norbp,isorb to match the total distributed scheme
     orbs%norbp = n_occ - orbsw%isorb
     if (orbsw%isorb + orbsw%norbp < n_occ ) orbs%norbp = orbsw%norbp
     if(orbsw%isorb > n_occ) orbs%norbp = 0
     orbs%isorb = orbsw%isorb
     call f_free_ptr(orbs%iokpt)
     orbs%iokpt = f_malloc_ptr(orbs%norbp,id='orbs%iokpt')
     orbs%iokpt=1
     if(orbs%norbp > 0) then
        nullify(orbs%eval)
        orbs%eval = f_malloc_ptr(orbs%norb*orbs%nkpts,id='orbs%eval')
        filename=trim(input%dir_output) // 'wavefunction'
        call readmywaves(iproc,filename,iformat,orbs,Glr%d%n1,Glr%d%n2,Glr%d%n3,&
             & input%hx,input%hy,input%hz,atoms,rxyz_old,atoms%astruct%rxyz,  & 
             Glr%wfd,psi(1,1))
        call f_free_ptr(orbs%eval)
     end if

     ! For the non-occupied orbitals, need to change norbp,isorb
     orbsv%norbp = orbsw%isorb + orbsw%norbp - n_occ
     if (orbsw%isorb + orbsw%norbp < n_occ ) orbsv%norbp = 0
     if (orbsw%isorb > n_occ) orbsv%norbp = orbsw%norbp
     orbsv%isorb = 0
     if(orbsw%isorb >= n_occ) orbsv%isorb = orbsw%isorb - n_occ
     if(associated(orbsv%iokpt)) then
        call f_free_ptr(orbsv%iokpt)
     end if
     orbsv%iokpt = f_malloc_ptr(orbsv%norbp,id='orbsv%iokpt')
     orbsv%iokpt=1

     ! read unoccupied wavefunctions
     if(orbsv%norbp > 0) then
        filename=trim(input%dir_output) // 'virtuals'
        if(associated(orbsv%eval)) nullify(orbsv%eval)
        orbsv%eval = f_malloc_ptr(orbsv%norb*orbsv%nkpts,id='orbsv%eval')
        call readmywaves(iproc,filename,iformat,orbsv,Glr%d%n1,Glr%d%n2,Glr%d%n3,&
             & input%hx,input%hy,input%hz,atoms,rxyz_old,atoms%astruct%rxyz,  & 
             Glr%wfd,psi(1,1+orbs%norbp),virt_list)
        call f_free_ptr(orbsv%eval)
     end if


     call f_free_ptr(rxyz_old)
     call f_free_ptr(virt_list)


     call timing(iproc,'CrtProjectors ','OF')

     if (iproc==0) then
        call yaml_map('Reading the wavefunctions',.true.)
        write(*,*) '!===================================!'
        !write(*,*) '!===================================!'
        !write(*,*) '!  Reading the wavefunctions : DONE !'
        !write(*,*) '!===================================!'
     end if

     !#########################################################
     ! Construct the Wannier functions
     !#########################################################

     if (iproc == 0) then
        call yaml_open_sequence('Constructing WFs')
        !write(*,*) '!==================================!'
        !write(*,*) '!     Constructing WFs :           !'
        !write(*,*) '!==================================!'
     end if


     allocate(wann(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f),stat=i_stat)
     call memocc(i_stat,wann,'wann',subname) 
     allocate(wannr(Glr%d%n1i*Glr%d%n2i*Glr%d%n3i),stat=i_stat)
     call memocc(i_stat,wannr,'wannr',subname) 
     call initialize_work_arrays_sumrho(Glr,w)


     ! Separate plotwann
!     allocate(plotwann_par(0:nproc-1),stat=i_stat)
!     call memocc(i_stat,plotwann_par,'plotwann_par',subname)
!     allocate(isplotwann_par(0:nproc-1),stat=i_stat)
!     call memocc(i_stat,isplotwann_par,'isplotwann_par',subname)
!     call parallel_repartition_with_kpoints(nproc,1,plotwann,plotwann_par)
!     ntot=0
!     do jproc=0,nproc-1
!        isplotwann_par(jproc)=ntot
!        ntot=ntot+plotwann_par(jproc)
!     end do
!     isplotwann = isplotwann_par(iproc) 
!     plotwannp = plotwann_par(iproc)
!     i_all = -product(shape(plotwann_par))*kind(plotwann_par)
!     deallocate(plotwann_par,stat=i_stat)
!     call memocc(i_stat,i_all,'plotwann_par',subname)
!     i_all = -product(shape(isplotwann_par))*kind(isplotwann_par)
!     deallocate(isplotwann_par,stat=i_stat)
!     call memocc(i_stat,i_all,'isplotwann_par',subname)


     ! Now construct the WFs
     ifile = 12 + iproc

     do iiwann = 1, nwannCon
        iwann = ConstList(iiwann)
        if (iwann == 0) then
           call yaml_warning('this should not happen')
           stop
        end if
        call to_zero(Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f, wann(1))
        do iband = 1, orbsw%norbp
           do i = 1, Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f
              wann(i) = wann(i) + umn(iwann,iband+orbsw%isorb) * psi(i,iband)
           end do
        end do

        ! Construction of the Wannier function.
        call mpiallred(wann(1),Glr%wfd%nvctr_c+7*Glr%wfd%nvctr_f,MPI_SUM)

        if(iproc == 0) then
           write(num,'(i4.4)') iwann
           if(outputype == 'cube') then
              !Put it in interpolating scaling functions
              call daub_to_isf(Glr,w,wann(1),wannr)
              call write_wannier_cube(ifile,trim(seedname)//'_'//num//'.cube',atoms,Glr,input,atoms%astruct%rxyz,wannr)
           else if(trim(outputype)=='bin') then
              call f_free_ptr(orbsw%iokpt)
              orbsw%iokpt = f_malloc_ptr(nwannCon,id='orbsw%iokpt')
              orbsw%iokpt=1
              if(hamilana .and. linear) then
                 call yaml_sequence(advance='no')
                 call yaml_open_map('orbs',flow=.true.)
                 call yaml_map('iiwann',iiwann)
                 call yaml_map('iwann',iwann)
                 call yaml_map('iokpt', (/ orbsw%iokpt(iiwann),orbsw%iokpt(iwann) /))
                 call yaml_close_map()
                 !write(*,*) 'iokpt',iiwann,orbsw%iokpt(iiwann),iwann,orbsw%iokpt(iwann)
                 call open_filename_of_iorb(ifile,.not.outformat,'minBasis',orbsw,iiwann,1,iwann_out)
                 ldim = Lzd%Llr(iiwann)%wfd%nvctr_c+7*Lzd%Llr(iiwann)%wfd%nvctr_f
                 gdim = Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f
                 allocate(lwann(ldim),stat=i_stat)
                 call memocc(i_stat,lwann,'lwann',subname)
                !DEBUG
                 !call plot_wf(trim(seedname)//'_'//num,1,atoms,1.0d0,Lzd%Glr,Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),rxyz,wann)
                 !call wpdot_wrap(1,Lzd%Glr%wfd%nvctr_c,Lzd%Glr%wfd%nvctr_f,Lzd%Glr%wfd%nseg_c,Lzd%Glr%wfd%nseg_f,&
                 !                 Lzd%Glr%wfd%keyvglob,Lzd%Glr%wfd%keyglob,wann,&
                 !                 Lzd%Glr%wfd%nvctr_c,Lzd%Glr%wfd%nvctr_f,Lzd%Glr%wfd%nseg_c,Lzd%Glr%wfd%nseg_f,&
                 !                 Lzd%Glr%wfd%keyvglob,Lzd%Glr%wfd%keyglob,wann,Gnorm)
                !END DEBUG
                 call psi_to_locreg2(iproc, ldim, gdim, Lzd%Llr(iiwann), Lzd%Glr, wann, lwann)
                !DEBUG
                 !call wpdot_wrap(1,Lzd%Llr(iiwann)%wfd%nvctr_c,Lzd%Llr(iiwann)%wfd%nvctr_f,Lzd%Llr(iiwann)%wfd%nseg_c,&
                 !                 Lzd%Llr(iiwann)%wfd%nseg_f,Lzd%Llr(iiwann)%wfd%keyvglob,Lzd%Llr(iiwann)%wfd%keyglob,lwann,&
                 !                 Lzd%Llr(iiwann)%wfd%nvctr_c,Lzd%Llr(iiwann)%wfd%nvctr_f,Lzd%Llr(iiwann)%wfd%nseg_c,&
                 !                 Lzd%Llr(iwann)%wfd%nseg_f,Lzd%Llr(iiwann)%wfd%keyvglob,Lzd%Llr(iiwann)%wfd%keyglob,lwann,&
                 !                 Lnorm)
                 !write(*,*) 'Norm of wann function',iwann, 'is:',Gnorm, 'while the cutting yields:',Lnorm
                 !call to_zero(gdim,wann(1))
                 !call Lpsi_to_global2(iproc, ldim, gdim, 1, 1, 1, Lzd%Glr,Lzd%Llr(iiwann), lwann(1), wann(1))
                 !Put it in interpolating scaling functions
                 !call daub_to_isf(Lzd%Glr,w,wann(1),wannr)
                 !call write_wannier_cube(ifile,trim(seedname)//'_test_'//num//'.cube',atoms,Glr,input,rxyz,wannr)
                 !stop 
                !END DEBUG
                 call writeonewave_linear(ifile,outformat,iiwann,Glr%d%n1,Glr%d%n2,Glr%d%n3,input%hx,input%hy,input%hz, &
                   (/cxyz(1,iwann),cxyz(2,iwann),cxyz(3,iwann) /),locrad(iwann),4,0.0d0,atoms%astruct%nat,atoms%astruct%rxyz,  & 
                   Lzd%Llr(iiwann)%wfd%nseg_c,Lzd%Llr(iiwann)%wfd%nvctr_c,Lzd%Llr(iiwann)%wfd%keyglob(1,1),&
                   Lzd%Llr(iiwann)%wfd%keyvglob(1),Lzd%Llr(iiwann)%wfd%nseg_f,Lzd%Llr(iiwann)%wfd%nvctr_f,&
                   Lzd%Llr(iiwann)%wfd%keyglob(1,Lzd%Llr(iiwann)%wfd%nseg_c+1),&
                   Lzd%Llr(iiwann)%wfd%keyvglob(Lzd%Llr(iiwann)%wfd%nseg_c+1), & 
                   lwann(1),lwann(Lzd%Llr(iiwann)%wfd%nvctr_c+1), ham(1,iwann,iwann))
                  i_all = -product(shape(lwann))*kind(lwann)
                  deallocate(lwann,stat=i_stat)
                  call memocc(i_stat,i_all,'lwann',subname)
              else if(hamilana) then
                ! open(ifile, file=trim(seedname)//'_'//num//'.bin', status='unknown',form='formatted')
                 call open_filename_of_iorb(ifile,.not.outformat,trim(seedname),orbsw,iwann,1,iwann_out)
                 call writeonewave(ifile,outformat,iiwann,Glr%d%n1,Glr%d%n2,Glr%d%n3,input%hx,input%hy,input%hz,&
                   atoms%astruct%nat,atoms%astruct%rxyz,  & 
                   Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keygloc(1,1),Glr%wfd%keyvloc(1),  & 
                   Glr%wfd%nseg_f,Glr%wfd%nvctr_f,Glr%wfd%keygloc(1,Glr%wfd%nseg_c+1),Glr%wfd%keyvloc(Glr%wfd%nseg_c+1), & 
                   wann(1),wann(Glr%wfd%nvctr_c+1), ham(1,iwann,iwann))
              else
                 call writeonewave(ifile,outformat,iiwann,Glr%d%n1,Glr%d%n2,Glr%d%n3,input%hx,input%hy,input%hz,&
                   atoms%astruct%nat,atoms%astruct%rxyz,  & 
                   Glr%wfd%nseg_c,Glr%wfd%nvctr_c,Glr%wfd%keygloc(1,1),Glr%wfd%keyvloc(1),  & 
                   Glr%wfd%nseg_f,Glr%wfd%nvctr_f,Glr%wfd%keygloc(1,Glr%wfd%nseg_c+1),Glr%wfd%keyvloc(Glr%wfd%nseg_c+1), & 
                   wann(1),wann(Glr%wfd%nvctr_c+1), -0.5d0)
              end if
              close(ifile)
           else
              stop 'ETSF not implemented yet'                
              ! should be write_wave_etsf  (only one orbital)
              call write_waves_etsf(iproc,trim(seedname)//'_'//num//'.etsf',orbs,Glr%d%n1,Glr%d%n2,Glr%d%n3,&
                   input%hx,input%hy,input%hz,atoms,atoms%astruct%rxyz,Glr%wfd,wann)
           end if
        end if
     end do  ! closing loop on iwann
!stop
     !Now if linear, write the coefficients
     if(linear .and. iproc == 0)then
        ! First construct the appropriate coefficient matrix
        allocate(umnt(nwannCon,nband),stat=i_stat)
        call memocc(i_stat,umnt,'umnt',subname)
        do iiwann = 1, nwannCon
           iwann = ConstList(iiwann)
           if(iwann == 0) then
              call yaml_warning('Umnt construction: this should not happen')
              stop
           end if
           do iiband = 1, nbandCon
              iband = bandlist(iiband)
              umnt(iiwann,iiband) = umn(iwann,iband) 
           end do
        end do
        ! write the coefficients to file
        if(outformat) then
           open(unit=99, file='minBasis'//'_coeff.bin', status='unknown',form='formatted')
        else 
           open(unit=99, file='minBasis'//'_coeff.bin', status='unknown',form='unformatted')
        end if
        call writeLinearCoefficients(99,outformat,Glr%d%n1,Glr%d%n2,Glr%d%n3,&
             & input%hx,input%hy,input%hz,atoms%astruct%nat,atoms%astruct%rxyz,&
             nbandCon,nwannCon,Glr%wfd%nvctr_c,Glr%wfd%nvctr_f,umnt)
        close(unit=99)
     end if

     if (iproc==0) then
        call yaml_close_sequence()
        !write(*,*) '!==================================!'
        !write(*,*) '!     Constructing WFs : DONE      !'
        !write(*,*) '!==================================!'
     end if


!DEBUG WRITE AND READ OF LINEAR WAVEFUNCTIONS
! Fake TMB orbs
!  call orbitals_descriptors(iproc,nproc,nwannCon,nwannCon,0, &
!       & orbsw%nspin,orbsw%nspinor,orbsw%nkpts,orbsw%kpts,orbsw%kwgts,wannorbs,.false.)
!  call orbitals_communicators(iproc,nproc,Glr,wannorbs,wanncomms)
!
! First must initialize Lzd (must contain Glr and hgrids)
!  call nullify_local_zone_descriptors(Lzd)
!  call copy_locreg_descriptors(Glr, Lzd%Glr, subname)
!  lzd%hgrids(1)=input%hx
!  lzd%hgrids(2)=input%hy
!  lzd%hgrids(3)=input%hz
!
! Then allocate some stuff
!  allocate(rxyz_old(3,atoms%astruct%nat),stat=i_stat)
!  call memocc(i_stat,rxyz_old,'rxyz_old',subname)
!
!  call initialize_linear_from_file(iproc,nproc,'minBasis',WF_FORMAT_BINARY,Lzd,wannorbs,atoms,rxyz)
!
!
!  do ilr=1,Lzd%nlr
!     if (iproc == 0) then
!        write(*,*)'Description of zone:',ilr
!        write(*,*)'ns:',Lzd%Llr(ilr)%ns1,Lzd%Llr(ilr)%ns2,Lzd%Llr(ilr)%ns3
!        write(*,*)'ne:',Lzd%Llr(ilr)%ns1+Lzd%Llr(ilr)%d%n1,Lzd%Llr(ilr)%ns2+Lzd%Llr(ilr)%d%n2,Lzd%Llr(ilr)%ns3+Lzd%Llr(ilr)%d%n3
!        write(*,*)'n:',Lzd%Llr(ilr)%d%n1,Lzd%Llr(ilr)%d%n2,Lzd%Llr(ilr)%d%n3
!        write(*,*)'nfl:',Lzd%Llr(ilr)%d%nfl1,Lzd%Llr(ilr)%d%nfl2,Lzd%Llr(ilr)%d%nfl3
!        write(*,*)'nfu:',Lzd%Llr(ilr)%d%nfu1,Lzd%Llr(ilr)%d%nfu2,Lzd%Llr(ilr)%d%nfu3
!        write(*,*)'ni:',Lzd%Llr(ilr)%d%n1i,Lzd%Llr(ilr)%d%n2i,Lzd%Llr(ilr)%d%n3i
!        write(*,*)'outofzone',ilr,':',Lzd%Llr(ilr)%outofzone(:)
!        write(*,*)'center: ',(Lzd%Llr(ilr)%locregCenter(i),i=1,3)
!        write(*,*)'locrad: ',Lzd%Llr(ilr)%locrad
!        write(*,*)'wfd dimensions: ',Lzd%Llr(ilr)%wfd%nseg_c, Lzd%Llr(ilr)%wfd%nseg_f, Lzd%Llr(ilr)%wfd%nvctr_c,&
!                   Lzd%Llr(ilr)%wfd%nvctr_f
!     end if
!  end do
!
!  call wavefunction_dimension(Lzd,wannorbs)
!  !npsidim = 0
!  !do iorb=1,orbs%norbp
!  !   ilr = orbs%inwhichlocreg(iorb+orbs%isorb)
!  !   npsidim = nspidim + Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f
!  !end do
!
!  allocate(psi2(wannorbs%npsidim_orbs),stat=i_stat)
!  call memocc(i_stat,psi,'psi',subname)
!  allocate(wannorbs%eval(wannorbs%norb),stat=i_stat)
!  call memocc(i_stat,wannorbs%eval,'wannorbs%eval',subname)
!  allocate(coeff(nbandCon,wannorbs%norb),stat=i_stat)
!  call memocc(i_stat,coeff,'coeff',subname)
! 
!  write(*,*) 'before reading the TMBs'
!  call readmywaves_linear(iproc,'minBasis',WF_FORMAT_BINARY,nbandCon,Lzd,wannorbs,atoms,rxyz_old,rxyz,  & 
!    psi2,coeff)
!  write(*,*) 'after reading the TMBs', sum(psi2)
!
!  write(*,*) 'The coefficients are:'
!  do i_all = 1, nbandCon
!    do ilr=1,wannorbs%norb
!       write(*,*) 'coeff',i_all, ilr, coeff(ilr,i_all)
!    end do
!  end do
!  i_all = -product(shape(rxyz_old))*kind(rxyz_old)
!  deallocate(rxyz_old,stat=i_stat)
!  call memocc(i_stat,i_all,'rxyz_old',subname)
!
!  indL = 1
!  do iwann = 1, wannorbs%norb
!     write(num,'(i4.4)') iwann
!     ilr = wannorbs%inwhichlocreg(iwann)
!     ldim = Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(iiwann)%wfd%nvctr_f
!     gdim = Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f
!     !allocate(wann(gdim),stat=i_stat)
!     !call memocc(i_stat,lwann,'lwann',subname)
!     call to_zero(gdim,wann(1))
!     call Lpsi_to_global2(iproc, ldim, gdim, 1, 1, 1, Lzd%Glr, Lzd%Llr(ilr),&
!          psi2(indL), wann(1))
!     indL = indL + ldim
!     !Put it in interpolating scaling functions
!     call daub_to_isf(Lzd%Glr,w,wann(1),wannr)
!     call write_wannier_cube(ifile,trim(seedname)//'_'//num//'.cube',atoms,Glr,input,rxyz,wannr)
!     !i_all = -product(shape(wann))*kind(wann)
!     !deallocate(wann,stat=i_stat)
!     !call memocc(i_stat,i_all,'lwann',subname)
!  end do
!
!END DEBUG

     if(linear)then
        call deallocate_local_zone_descriptors(Lzd, subname)
        i_all = -product(shape(locrad))*kind(locrad)
        deallocate(locrad,stat=i_stat)
        call memocc(i_stat,i_all,'locrad',subname)
        i_all = -product(shape(cxyz))*kind(cxyz)
        deallocate(cxyz,stat=i_stat)
        call memocc(i_stat,i_all,'cxyz',subname)
        if(iproc == 0) then
           i_all = -product(shape(umnt))*kind(umnt)
           deallocate(umnt,stat=i_stat)
           call memocc(i_stat,i_all,'umnt',subname)
        end if
        i_all = -product(shape(calcbounds))*kind(calcbounds)
        deallocate(calcbounds,stat=i_stat)
        call memocc(i_stat,i_all,'calcbounds',subname)
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
  i_all = -product(shape(Constlist))*kind(Constlist)
  deallocate(Constlist,stat=i_stat)
  call memocc(i_stat,i_all,'Constlist',subname)
  i_all = -product(shape(bandlist))*kind(bandlist)
  deallocate(bandlist,stat=i_stat)
  call memocc(i_stat,i_all,'bandlist',subname)
  i_all = -product(shape(umn))*kind(umn)
  deallocate(umn,stat=i_stat)
  call memocc(i_stat,i_all,'umn',subname)
  if(bondana .or. hamilana) then
     i_all = -product(shape(buf))*kind(buf)
     deallocate(buf,stat=i_stat)
     call memocc(i_stat,i_all,'buf',subname)
  end if
  if(hamilana) then
     i_all=-product(shape(ham))*kind(ham)
     deallocate(ham,stat=i_stat)
     call memocc(i_stat,i_all,'ham',subname)
  end if
  i_all = -product(shape(wann_list))*kind(wann_list)
  deallocate(wann_list,stat=i_stat)
  call memocc(i_stat,i_all,'wann_list',subname)
  call deallocate_lr(Glr,subname)
  call deallocate_orbs(orbs,subname)
  !call deallocate_atoms_scf(atoms,subname)
  call deallocate_atoms_data(atoms)
!  call free_input_variables(input)
  call free_input_variables(input)
  call f_lib_finalize()
  !free all yaml_streams active
  call yaml_close_all_streams()

  !#########################################################
  ! Ending timing and MPI
  !#########################################################
  call f_timing_stop(mpi_comm=bigdft_mpi%mpi_comm,nproc=bigdft_mpi%nproc,gather_routine=gather_timings)

  call cpu_time(tcpu1)
  call system_clock(ncount1,ncount_rate,ncount_max)
  tel=real(ncount1-ncount0,kind=8)/real(ncount_rate,kind=8)
  if (iproc == 0) &
    call yaml_map('CPU time/ELAPSED time for root process ', (/ tel,real(tcpu1-tcpu0,kind=8) /),fmt='(f12.2)')
    !write( *,'(1x,a,1x,i4,2(1x,f12.2))') 'CPU time/ELAPSED time for root process ', iproc,tel,tcpu1-tcpu0
   !finalize memory counting
!   call memocc(0,0,'count','stop')
 
  call bigdft_finalize(ierr)
!
!  ! Barrier suggested by support for titane.ccc.cea.fr, before finalise.
!  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
! 
!  call MPI_FINALIZE(ierr)

end program Waco

subroutine Waco_input_variables(iproc,filename,nband,nwann,bondAna,Stereo,hamilAna,WannCon,filetype,nwannCon,refpos,units,&
           sprdfact,sprddiff,enediff,iformat,linear,nbandmB,sprdmult)
   use module_base
   use module_types
   use module_input
   implicit none
   integer, intent(in) :: iproc
   character(len=*), intent(in) :: filename
   logical, intent(out) :: bondAna, Stereo, hamilAna, WannCon, iformat, linear
   character(len=4), intent(out) :: filetype,units
   integer, intent(out) :: nband,nwann,nwannCon, nbandmB
   real(gp), intent(out) :: sprdfact,sprddiff,enediff,sprdmult
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
   call input_var(filetype,'cube')
   call input_var(iformat,'F', comment='! Wannier function construction, type of output file (cube, bin or etsf)')
   call input_var(linear,'F')
   call input_var(nbandmB,'1')
   call input_var(sprdmult,'3.0',comment='!Ouput as minimal basis, spread factor for minimal basis construction')
   call input_var(nwannCon,'1',comment='! number of Wannier to construct (if 0 do all) followed by Wannier list on next &
   &     line (optional)')
   
   !Make defaults
   call input_free()

end subroutine Waco_input_variables

subroutine read_input_waco(filename,nwannCon,Constlist,linear,nbandmB,bandlist)
  use module_base
  implicit none
  character(len=*), intent(in) :: filename
  logical, intent(in) :: linear
  integer, intent(in) :: nwannCon,nbandmB
  integer, dimension(nwannCon), intent(out) :: Constlist
  integer, dimension(nbandmB), intent(out) :: bandlist
 
  ! Local variable
  integer :: i,ierr

  open(unit=22, file=trim(filename),status='old')
  ! Skip first lines
  do i=1,7
     read(22,*)
  end do

  ! read the Constlist
  read(22,*,iostat=ierr) (Constlist(i),i=1,nwannCon)
  close(unit=22)
  if(ierr < 0) then  !reached the end of file and no Constlist, so generate the trivial one
    do i= 1, nwannCon
       Constlist(i) = i
    end do
  end if

  !read the bandlist
  if(linear)then
    read(22,*,iostat=ierr) (bandlist(i),i=1,nbandmB)
    if(ierr < 0) then  !reached the end of file and no bandlist, so generate the trivial one
      do i= 1, nbandmB
         bandlist(i) = i
      end do
    end if
  else
    bandlist = 0
  end if
  

end subroutine read_input_waco


!> This routine reads the list of virtual orbitals needed
subroutine read_inter_list(iproc,n_virt, virt_list)

   use yaml_output

   implicit none

   ! I/O variables
   integer, intent(in) :: n_virt,iproc
   integer, dimension(n_virt), intent(out) :: virt_list

   ! Local variables
   character(len=*), parameter :: filename='input.inter'
   integer :: i,j,ierr 

   open(unit=11, file=filename, status='old')

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
   close(unit=11)

   if (iproc==0) then
      call yaml_map('Reading virtual orbitals',.true.)
      !write(*,*) '!==================================!'
      !write(*,*) '!Reading virtual orbitals list done!'
      !write(*,*) '!==================================!'
   end if

END SUBROUTINE read_inter_list


!> This routine reads the first lines of a .inter file
subroutine read_inter_header(iproc,seedname, filetype,residentity,write_resid, n_occ, pre_check,&
           n_virt_tot, n_virt, w_unk, w_sph, w_ang, w_rad)

   use yaml_output
   implicit none

   ! I/O variables
   integer, intent(in) :: iproc
   character, intent(out) :: seedname*16, filetype*4
   integer, intent(out) :: n_occ, n_virt, n_virt_tot
   logical, intent(out) :: w_unk, w_sph, w_ang, w_rad, pre_check, residentity, write_resid

   ! Local variables
   character(len=*), parameter :: filename='input.inter'
   character :: char1*1, char2*1, char3*1, char4*1
   logical :: file_exist
   integer :: ierr

   ! Should check if it exists, if not, make a nice output message
   inquire(file=filename,exist=file_exist)
   if (.not. file_exist) then
      if(iproc == 0) then
         call yaml_warning('Input file, input.inter, not found !')
         call yaml_warning('CORRECTION: Create or give correct input.inter file.')
         !write(*,'(A)') 'ERROR : Input file, input.inter, not found !'
         !write(*,'(A)') 'CORRECTION: Create or give correct input.inter file.'
      end if
      call mpi_finalize(ierr)
      stop
   end if

   open(unit=11, FILE=filename, STATUS='OLD')

   if (iproc == 0) then
      call yaml_open_map('Reading input.inter header')
      !write(*,*) '!==================================!'
      !write(*,*) '!   Reading input.inter header :   !'
      !write(*,*) '!==================================!'
   end if

   w_unk=.false.
   w_sph=.false.
   w_ang=.false.
   w_rad=.false.

   ! First line
   read(11,*) seedname
   if (iproc == 0) call yaml_map('System studied', trim(seedname))
   !if (iproc == 0) write(*,*) 'System studied : ', trim(seedname)

   ! Second line
   read(11,*) filetype
   if (iproc == 0) call yaml_map('File type', trim(filetype))
   !if (iproc == 0) write(*,*) 'file type : ', filetype

   ! Third line
   read(11,*) char1, char2, n_occ
   if (iproc == 0) call yaml_map('Number of occupied orbitals', n_occ)
   !if (iproc == 0) write(*,'(A30,I4)') 'Number of occupied orbitals :', n_occ
   if(char1=='T') then
     residentity = .true.
     if (iproc == 0) call yaml_map('Use resolution of the identity to construct virtual states',.true.)
     !if (iproc == 0) write(*,*) 'Will use resolution of the identity to construct virtual states'
   else
     residentity = .false.
   end if
   if(residentity .and. char2=='T')then
     write_resid = .true.
     if (iproc == 0) call yaml_map('The constructed virtual states will be written to file',trim(filename))
     !if (iproc == 0) write(*,*) 'The constructed virtual states will be written to file.'
   else
     write_resid = .false.
   end if

   ! Fourth line
   read(11,*) char1, n_virt_tot, n_virt
   if (char1=='T' .and. .not. residentity) then
      pre_check=.true.
      !if (iproc == 0) write(*,*) 'Pre-check before calculating Amnk and Mmnk matrices'
      !if (iproc == 0) write(*,'(A38,I4)') 'Total number of unoccupied orbitals :', n_virt_tot
   else if(.not. residentity)then
      pre_check=.false.
      !if (iproc == 0) write(*,*) 'Calculation of Amnk and Mmnk matrices'
      !if (iproc == 0) write(*,'(A39,I4)') 'Number of chosen unoccupied orbitals :', n_virt
   else
      pre_check=.false.
      !if (iproc == 0) write(*,'(A38,I4)') 'Total number of unoccupied orbitals :', n_virt_tot
   end if
   if (iproc == 0) then
      call yaml_map('Pre-check before calculating Amnk and Mmnk matrices',pre_check)
      call yaml_map('Total number of unoccupied orbitals', n_virt_tot)
      call yaml_map('Number of chosen unoccupied orbitals :', n_virt)
   end if

   ! Fifth line
   read(11,*) char1, char2, char3, char4
   if (char1=='T') then
      w_unk=.true.
      if (iproc == 0) call yaml_comment('You want to write a UNKp.s file')
      !if (iproc == 0) write(*,*) 'You want to write a UNKp.s file'
   else if (char1 /= 'F') then
      if (iproc == 0) call yaml_warning('Wrong value for w_unk')
      !if (iproc == 0) write(*,*) 'Wrong value for w_unk'
      STOP
   end if
   if (char2=='T') then
      w_sph=.true.
      if (iproc == 0) call yaml_comment('You want to write .cube files for spherical harmonics')
      !if (iproc == 0) write(*,*) 'You want to write .cube files for spherical harmonics'
   else if (char2 .ne. 'F') then
      if (iproc == 0) call yaml_warning('Wrong value for w_sph')
      !if (iproc == 0) write(*,*) 'Wrong value for w_sph'
      STOP
   end if
   if (char3=='T') then
      w_ang=.true.
      if (iproc == 0) call yaml_comment('You want to write .cube files for angular parts of the spherical harmonics')
      !if (iproc == 0) write(*,*) 'You want to write .cube files for angular parts of the spherical harmonics'
   else if (char3 .ne. 'F') then
      if (iproc == 0) call yaml_warning('Wrong value for w_ang')
      !if (iproc == 0) write(*,*) 'Wrong value for w_ang'
      STOP
   end if
   if (char4=='T') then
      w_rad=.true.
      if (iproc == 0) call yaml_comment('You want to write .cube files for radial parts of the spherical harmonics')
      !if (iproc == 0) write(*,*) 'You want to write .cube files for radial parts of the spherical harmonics'
   else if (char4 .ne. 'F') then
      if (iproc == 0) call yaml_warning('Wrong value for w_rad')
      !if (iproc == 0) write(*,*) 'Wrong value for w_rad'
      STOP
   end if

   close(unit=11)

   if (iproc == 0) then
      call yaml_close_map() 
      !write(*,*) '!==================================!'
      !write(*,*) '! Reading input.inter header done  !'
      !write(*,*) '!==================================!'
   end if

end subroutine read_inter_header


subroutine read_umn(iproc,nwann,nband,seedname,umn)
   use module_types
   use yaml_output
   implicit none
   integer, intent(in) :: iproc
   integer, intent(in) :: nwann, nband
   character(len=*),intent(in) :: seedname
   real(gp),dimension(nwann,nband),intent(out) :: umn
   !Local variables
   logical :: file_exist
   integer :: ierr, nwann_umn,nband_umn,iwann,iband,int1,int2

   if (iproc == 0) then
      call yaml_open_map('Reading .umn')
      write(*,*) '!==================================!'
      !write(*,*) '!==================================!'
      !write(*,*) '!     Reading .umn :               !' 
      !write(*,*) '!==================================!'
   end if

   inquire(file=trim(seedname)//'.umn',exist=file_exist)
   if (.not. file_exist) then
      if (iproc==0) then
         call yaml_warning('Input file,' // trim(seedname) // '.umn, not found !')
         call yaml_warning('CORRECTION: Create or give correct input file.')
         !write(*,'(A,1x,A)') 'ERROR : Input file,',trim(seedname)//'.umn, not found !'
         !write(*,'(A)') 'CORRECTION: Create or give correct input file.'
      end if
      call mpi_finalize(ierr)
      stop
   end if

   open(unit=11, file=trim(seedname)//'.umn', status='OLD')
   read(11,*) nwann_umn, nband_umn

   if(nwann_umn .ne. nwann .or. nband_umn .ne. nband) then
     if(iproc == 0) then
       call yaml_warning('Number of wannier functions in umn,' // trim(yaml_toa(nwann_umn)) // &
          & 'not equal number of Wannier functions used:' // trim(yaml_toa(nwann)))
       call yaml_warning('Number of orbitals in umn,' // trim(yaml_toa(nband_umn)) // &
          & 'not equal number of orbitals used:' // trim(yaml_toa(nband)))
       !write(*,'(A,I4)') 'ERROR : number of wannier functions in umn,',nwann_umn
       !write(*,'(A,I4)') 'not equal number of Wannier functions used:',nwann
       !write(*,'(A,I4)') 'ERROR : number of orbitals in umn,',nband_umn
       !write(*,'(A,I4)') 'not equal number of orbitals used ',nband
     end if
     call mpi_finalize(ierr)
     stop 
   end if

   do iwann = 1,nwann
      do iband = 1,nband
         read(11,*) int1, int2, umn(iwann,iband)
      end do
   end do
   close(unit=11)

   if (iproc == 0) then
      call yaml_map('Number of Wannier functions',nwann_umn)
      call yaml_map('Number of orbitals used',nband_umn)
      call yaml_close_map()
      !write(*,*) 'Number of Wannier functions : ',nwann_umn
      !write(*,*) 'Number of orbitals used     : ',nband_umn
      !write(*,*) '!==================================!'
      !write(*,*) '!     Reading .umn : DONE          !' 
      !write(*,*) '!==================================!'
   end if
end subroutine read_umn


subroutine read_centers(iproc,nwann,plotwann,natom,seedname,wann_list,cxyz,rxyz_wann,readsprd,sprd)
   use module_types
   use yaml_output
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
         call yaml_warning('Input file,' // trim(seedname) //'_centres.xyz, not found!')
         call yaml_warning('CORRECTION: Create or give correct input file.')
         !write(*,'(A,1x,A)') 'ERROR : Input file,',trim(seedname)//'_centres.xyz, not found !'
         !write(*,'(A)') 'CORRECTION: Create or give correct input file.'
      end if
      call mpi_finalize(i)
      stop
   end if

   open(unit=11, file=trim(seedname)//'_centres.xyz', status='OLD')

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
      iiwann = 0
      do iwann = 1, nwann
         commented = .true.
         do i = 1, plotwann
            if(iwann == wann_list(i)) commented = .false.
         end do
         if(.not. commented) then
           iiwann = iiwann + 1
           read(11, *) char1, sprd(iiwann)
         else
           read(11, *) ! just skip line
         end if
      end do
   end if
   close(unit=11)
end subroutine read_centers


subroutine read_nrpts_hamiltonian(iproc,seedname,nrpts)
   use yaml_output
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
         call yaml_warning('Input file,' // trim(seedname) // '_hr.dat, not found!')
         call yaml_warning('CORRECTION: Create or give correct input file.')
         !write(*,'(A,1x,A)') 'ERROR : Input file,',trim(seedname)//'_hr.dat, not found !'
         !write(*,'(A)') 'CORRECTION: Create or give correct input file.'
      end if
      call mpi_finalize(ierr)
      stop
   end if
   
   !read the hamiltonian
   open(unit=11, file=trim(seedname)//'_hr.dat', status='OLD')
   
   !skip first lines
   read(11,*)           !comment line containing the date 
   read(11,*)           !this line contains the number of wannier functions (already known)
   read(11,*) nrpts     ! this line contains the number of Wigner-Seitz grid points
   close(unit=11)

end subroutine read_nrpts_hamiltonian


subroutine read_hamiltonian(iproc,nrpts,nwann,seedname,ham)
   use module_types
   use yaml_output
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
         call yaml_warning('Input file,' // trim(seedname)//'_hr.dat, not found !')
         call yaml_warning('CORRECTION: Create or give correct input file.')
         !write(*,'(A,1x,A)') 'ERROR : Input file,',trim(seedname)//'_hr.dat, not found !'
         !write(*,'(A)') 'CORRECTION: Create or give correct input file.'
      end if
      call mpi_finalize(ierr)
      stop
   end if
   
   !read the hamiltonian
   open(unit=11, file=trim(seedname)//'_hr.dat', status='OLD')
   
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
   close(unit=11)
   
end subroutine read_hamiltonian

subroutine write_wannier_cube(ifile,filename,atoms,Glr,input,rxyz,wannr)
   use module_types
   implicit none
   character(len=*), intent(in) :: filename
   integer, intent(in) :: ifile
   type(atoms_data),intent(in) :: atoms
   type(locreg_descriptors), intent(in) :: Glr
   type(input_variables),intent(in) :: input
   real(gp),dimension(3,atoms%astruct%nat),intent(in) :: rxyz
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
   
   open(unit=ifile, file=filename, status='unknown')
   write(ifile,*) ' CUBE file for ISF field'
   write(ifile,*) ' Case for'
   write(ifile,'(I4,1X,F12.6,2(1X,F12.6))') atoms%astruct%nat, real(0.d0), real(0.d0), real(0.d0)
   write(ifile,'(I4,1X,F12.6,2(1X,F12.6))') Glr%d%n1i-(nbl1+nbr1), 0.5_dp*input%hx, real(0.d0),  real(0.d0)
   write(ifile,'(I4,1X,F12.6,2(1X,F12.6))') Glr%d%n2i-(nbl2+nbr2), real(0.d0),  0.5_dp*input%hy, real(0.d0)
   write(ifile,'(I4,1X,F12.6,2(1X,F12.6))') Glr%d%n3i-(nbl3+nbr3), real(0.d0),  real(0.d0),  0.5_dp*input%hz
   do i=1, atoms%astruct%nat
      write(ifile,'(I4,1X,F12.6,3(1X,F12.6))') atoms%nzatom(atoms%astruct%iatype(i)), real(0.d0), (real(rxyz(j,i)), j=1,3)
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
   close(unit=ifile)

end subroutine write_wannier_cube


subroutine scalar_kmeans_diffIG(iproc,nIG,crit,nel,vect,string,nbuf,buf)
  use BigDFT_API
  use module_interfaces
  use yaml_output
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
     call yaml_map('Number of iterations to reach the convergence',iter)
     call yaml_map( 'Number of elments for the clustering of ' // trim(string),nbuf)
     !write(*,'(A,1x,i4,1x,A)') 'Convergence reached in',iter,'iterations.'
     !write(*,'(A,A,A,1x,i4,1x,A)') 'The ',trim(string),' can be clustered in',nbuf,'elements:'
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

subroutine stereographic_projection(mode,natom, rxyz, refpos, CM, rad, proj, normal, dcp)
   use BigDFT_API
   use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
   use module_interfaces
   implicit none
   integer, intent(in) :: mode        ! 0= atomic projection, 1=wannier projection
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
!   open(unit=22,file='pos_sph.xyz', status='unknown')
!   write(22,'(I4)') natom+1
!   write(22,*) !skip this line
!   write(22,'(A,3(2x,E14.6))'), 'X',(refpos(j),j=1,3)
!   do iat = 1, natom
!      write(22,'(A,3(2x,E14.6))'),'B',pos(iat,1),pos(iat,2),pos(iat,3) 
!   end do
!   close(unit=22)
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
        if(dcp .ne. 0 .and. mode == 0) stop 'This should not happen'
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

!   do iat = 1, natom
!      do j= 1, 3
!         ! Now shift the positions to be contained in the first quadrant
!         proj2(iat,j) = proj(iat,j)! - minval(proj(:,j))
!      end do    
!   end do


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
   real(gp), dimension(atoms%astruct%nat,3), intent(in) :: proj  ! atom positions in the projection
   !Local variable
   integer :: i, iat

   ! Now open file that will contain the stereographic projection
   open(unitnumb,file=trim(filename), status='unknown')
   if(NeglectPoint .ne. 0) then
      write(unitnumb,*) atoms%astruct%nat-1
   else
      write(unitnumb,*) atoms%astruct%nat
   end if
   write(unitnumb,*)

   do iat = 1, atoms%astruct%nat
      ! Now print the information
      if(iat == NeglectPoint) then
         write(unitnumb,'(A,3(2x,E14.6))')'#   '//atoms%astruct%atomnames(atoms%astruct%iatype(iat)),(proj(iat,i), i=1,3)
      else
         write(unitnumb,'(A,3(2x,E14.6))')atoms%astruct%atomnames(atoms%astruct%iatype(iat)),(proj(iat,i), i=1,3)
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

subroutine build_stereographic_graph_facets(natoms,nsurf, mcenters,maxbond,rxyz,ncenters,Zatoms,nfacets,facets,vertex)
   use module_interfaces
   use module_types
   implicit none
   integer, intent(in) :: natoms, nsurf,mcenters
   real(gp), intent(in) :: maxbond
   real(gp), dimension(3,natoms), intent(in) :: rxyz
   integer, dimension(nsurf), intent(in) :: ncenters
   integer, dimension(natoms,nsurf),intent(in) :: Zatoms
   integer, dimension(nsurf), intent(out) :: nfacets
   integer, dimension(nsurf,mcenters*(mcenters - 1)/2, 3),intent(out) :: facets
   integer, dimension(nsurf,mcenters*(mcenters - 1)/2, 3),intent(out) :: vertex    !exactly like facets, but in reference to only the surface
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

subroutine output_stereographic_graph(natoms,mcenters,proj,projC,nsurf,ncenters,Zatoms,npoly,poly,vertex,normal,NeglectPoint) 
   use module_types
   implicit none
   integer, intent(in) :: natoms, mcenters, nsurf, NeglectPoint
   real(gp), dimension(natoms,3), intent(in) :: proj                   ! atom position in the projection
   real(gp), dimension(nsurf, 3), intent(in) :: projC                  ! Wannier centers in the projection
   integer, dimension(nsurf), intent(in) :: ncenters, npoly
   integer, dimension(natoms,nsurf),intent(in) :: Zatoms               ! indexes of all the atoms spanned by the Wannier function
   integer, dimension(nsurf,mcenters*(mcenters-1)/2,3),intent(in) :: poly
   integer, dimension(nsurf,mcenters*(mcenters-1)/2,3),intent(in) :: vertex
   real(gp), dimension(3), intent(in) :: normal
   ! Local variables
   integer :: isurf,ipts,i, nsurftot,npts,ndecimal,ncent, num_poly, num_poly_tot
   character(len=14) :: surfname
   character(len=20) :: forma
   logical :: condition

   open(unit=23,file='proj.surf', status='unknown')

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
   close(unit=23)

end subroutine output_stereographic_graph

subroutine wannier_dos(filename,nrpts,nwann,ntypes,wtypes,diag)
use module_types
implicit none
character(len=*), intent(in) :: filename
integer, intent(in) :: nrpts                           !< Number of unit cells (r-points)
integer, intent(in) :: nwann                           !< Number of Wannier functions 
integer, intent(in) :: ntypes                          !< Number of types of Wannier functions
integer, dimension(nrpts,nwann), intent(in) :: wtypes  !< Types of the Wannier functions
real(gp), dimension(nrpts,nwann),intent(in) :: diag    !< Diagonal elements of the Hamiltonian in Wannier basis
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
         ityp = wtypes(irpts,iwann) + 1
         dos(ityp,ipt) = dos(ityp,ipt) + gauss(ipt)
      end do
   end do
end do
write(numb,'(I3)') ntypes + 1
forma = '(E14.6,2x,'//numb//'E14.6)'

open(unit=22, file=trim(filename), status='unknown')
write(22,'(A)') '# Wannier DOS file'
write(22,'(A,2x,i4)') '# Number of Wannier functions:',nwann
write(22,'(A,2x,i4)') '# Number of real space points:',nrpts
write(22,'(A,2x,i4)') '# Number of energy points:',ndos
write(22,'(A,2x,E14.6)') '# Width of the Gaussians:', width
do ipt = 1, ndos
   ener = emin + real(ipt-51,kind=8)*inc
   write(22,forma) ener, (dos(i,ipt), i=1,ntypes+1)
end do
close(unit=22)

open(unit=22, file='integrated_'//trim(filename), status='unknown')
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
close(unit=22)


end subroutine wannier_dos

subroutine wannier_projected_dos(filename,nrpts,nwann,norb,umn,ntypes,wtypes,eigen)
use module_types
implicit none
character(len=*), intent(in) :: filename
integer, intent(in) :: nrpts                           !< Number of unit cells (r-points)
integer, intent(in) :: nwann                           !< Number of Wannier functions 
integer, intent(in) :: ntypes                          !< Number of types of Wannier functions
integer, intent(in) :: norb                            !< Number of orbitals
integer, dimension(nrpts,nwann), intent(in) :: wtypes  !< Types of the Wannier functions
real(gp), dimension(nwann,norb), intent(in) :: umn     !< Wannier transformation matrix
real(gp), dimension(norb),intent(in) :: eigen          !< Diagonal elements of the Hamiltonian in Wannier basis
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
            ityp = iwann+1 !wtypes(irpts,iwann) + 1
            dos(ityp,ipt) = dos(ityp,ipt) + gauss(ipt)
         end do
      end do
   end do
end do
write(numb,'(I3)') nwann +1
forma = '(E14.6,2x,'//numb//'E14.6)'

open(unit=22, file=trim(filename), status='unknown')
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
close(unit=22)

open(unit=22, file='integrated_'//trim(filename), status='unknown')
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
close(unit=22)


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

open(unit=22,file=filename,status='old')
do ikpt = 1, nkpt
   do iorb = 1, norb
      read(22,*) iiorb, iikpt, eigen(ikpt,iorb)
   end do
end do
close(unit=22)

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

open(unit=22,file=trim(filename)//'.amn',status='old')
read(22,*) ! skip first line which is a comment
read(22,*) nband, nkpt, nproj
close(unit=22)
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
integer :: int2, int3, int4, nk, np, nb
integer :: ierr
real(gp) :: r1

inquire(file=trim(filename)//'.amn',exist=file_exist)
if (.not. file_exist) then
   write(*,'(A)') 'ERROR : Input file,',trim(filename)//'.amn',', not found !'
   write(*,'(A)') 'CORRECTION: Create or give correct input.inter file.'
   call mpi_finalize(ierr)
   stop
end if

open(unit=22,file=trim(filename)//'.amn',status='old')
read(22,*) ! skip first line which is a comment
read(22,*) !skip this line: nband, nkpt, nproj
do nk=1, nkpt 
   do np=1, nproj
      do nb=1, nband
         read(22,*) int2, int3, int4, amn(nb,np), r1
      end do
   end do
end do
close(unit=22)

end subroutine read_amn


!> This routine reads an .nnkp file and returns the types of projectors
subroutine read_proj(seedname, n_kpts, n_proj, l, mr)
   implicit none
   ! I/O variables
   character(len=16),intent(in) :: seedname
   integer, intent(in) :: n_kpts, n_proj
   integer, dimension(n_proj), intent(out) :: l, mr
   ! Local variables
   integer :: i, j
   character(len=16) :: char1, char2, char3, char4
   logical :: calc_only_A
   real :: real_latt(3,3), recip_latt(3,3)
   real, dimension(n_kpts,3) :: kpts
   real(kind=8), dimension(n_proj,3) :: ctr_proj, x_proj, z_proj
   integer, dimension(n_proj) :: rvalue
   real, dimension(n_proj) :: zona


   open(unit=11, FILE=trim(seedname)//'.nnkp', STATUS='OLD')
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

   close(unit=11)

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
   character(len=2) :: num
   character(len=27):: forma
   character(len=10), dimension(nproj) :: label
   integer :: np, np2, iwann,iiwann, ntype, ii, i_stat, i_all
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
!   if (iproc == 0) then
     write(*,*) 'Analysis of the symmetry types of the Wannier functions'
     write(num,'(I2)') ntype
     forma = '(15x,'//trim(num)//'(A,6x))'
     write(*,trim(forma))(Wplabel(ii),ii=1,ntype)
     forma = '(I3,2x,I3,2x,'//trim(num)//'(E14.6,2x))'
     do iwann = 1, plotwann
        iiwann = wann_list(iwann)
        write(*,trim(forma)) iiwann, iwann, (Wpweight(iiwann,ii)/norm(iiwann), ii=1,ntype)
     end do
!   end if

    i_all = -product(shape(norm))*kind(norm)
    deallocate(norm,stat=i_stat)
    call memocc(i_stat,i_all,'norm',subname)
    i_all = -product(shape(Wpweight))*kind(Wpweight)
    deallocate(Wpweight,stat=i_stat)
    call memocc(i_stat,i_all,'Wpweight',subname)
    i_all = -product(shape(Wplabel))*kind(Wplabel)
    deallocate(Wplabel,stat=i_stat)
    call memocc(i_stat,i_all,'Wplabel',subname)

end subroutine character_list


subroutine read_spread_file(iproc,seedname,nwann,plotwann,wann_list)
implicit none
character(len=16),intent(in) :: seedname
integer, intent(in) :: iproc, nwann
integer, intent(out) :: plotwann
integer, dimension(nwann),intent(out) :: wann_list
!Local variables
character(len=8) :: string
logical :: file_exist
integer :: i,j

   if (iproc == 0) then
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

   !Initialize wann_list
   wann_list = 0

   open(unit=11, file=trim(seedname)//'.dat', status='OLD')

   ! Check if the line is commented
   plotwann = 0
   do i = 1 ,nwann
      read(11,*) string
      if(VERIFY(string, ' #', .false.) == 0) cycle
      plotwann = plotwann +1
      wann_list(plotwann) = i
   end do

   close(unit=11)

   if(iproc == 0) then
      write(*,'(A)') 'Occupied/Plotted Wannier functions are:'
      do i=1, plotwann, 6
         write(*,'(6(2x,I4))') (wann_list(j), j=i,min(i+5,size(wann_list)))
      end do
   end if

   if (iproc == 0) then
      write(*,*) '!==================================!'
      write(*,*) '!     Reading .dat : DONE          !'
      write(*,*) '!==================================!'
   end if
END SUBROUTINE read_spread_file


subroutine get_mindist(geocode,rxyz,cxyz,box,distw)
   use module_types
   implicit none
   character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
   real(gp),dimension(3), intent(in) :: rxyz, cxyz, box
   real(gp),intent(out) :: distw
   !Local variables
   integer :: i
   real(gp) :: dist1, dist2,dist3
   real(gp),dimension(3) :: distp

   if(geocode == 'F') then
     do i=1,3
        distp(i) = (rxyz(i)-cxyz(i))**2
     end do
   else if(geocode == 'S') then
      do i=1,3
         dist1 = (rxyz(i)-cxyz(i))**2
         if(i/=2) then
            dist2 = (rxyz(i)-cxyz(i)-box(i))**2
            dist3 = (rxyz(i)-cxyz(i)+box(i))**2
         else
            dist2 = box(i)
            dist3 = box(i)
         end if
         distp(i) = min(dist1,dist2,dist3)
      end do
   else if(geocode == 'P') then
      do i=1,3
         dist1 = (rxyz(i)-cxyz(i))**2
         dist2 = (rxyz(i)-cxyz(i)-box(i))**2
         dist3 = (rxyz(i)-cxyz(i)+box(i))**2
         distp(i) = min(dist1,dist2,dist3)
      end do
   end if

   distw = sqrt(distp(1) + distp(2) + distp(3))

END SUBROUTINE get_mindist 

!subroutine writeonewave_linear(unitwf,useFormattedOutput,iorb,n1,n2,n3,hx,hy,hz,locregCenter,&
!     locrad,confPotOrder,confPotprefac,nat,rxyz, nseg_c,nvctr_c,keyg_c,keyv_c,  &
!     nseg_f,nvctr_f,keyg_f,keyv_f, &
!     psi_c,psi_f,eval)
!  use module_base
!  implicit none
!  logical, intent(in) :: useFormattedOutput
!  integer, intent(in) :: unitwf,iorb,n1,n2,n3,nat,nseg_c,nvctr_c,nseg_f,nvctr_f,confPotOrder
!  real(gp), intent(in) :: hx,hy,hz,locrad,confPotprefac
!  real(wp), intent(in) :: eval
!  integer, dimension(nseg_c), intent(in) :: keyv_c
!  integer, dimension(nseg_f), intent(in) :: keyv_f
!  integer, dimension(2,nseg_c), intent(in) :: keyg_c
!  integer, dimension(2,nseg_f), intent(in) :: keyg_f
!  real(wp), dimension(nvctr_c), intent(in) :: psi_c
!  real(wp), dimension(7,nvctr_f), intent(in) :: psi_f
!  real(gp), dimension(3,nat), intent(in) :: rxyz
!  real(gp), dimension(3), intent(in) :: locregCenter
!  !local variables
!  integer :: iat,jj,j0,j1,ii,i0,i1,i2,i3,i,iseg,j
!  real(wp) :: tt,t1,t2,t3,t4,t5,t6,t7
!
!  if (useFormattedOutput) then
!     write(unitwf,*) iorb,eval
!     write(unitwf,*) hx,hy,hz
!     write(unitwf,*) n1,n2,n3
!     write(unitwf,*) locregCenter(1),locregCenter(2),locregCenter(3),locrad,confPotOrder, confPotprefac
!     write(unitwf,*) nat
!     do iat=1,nat
!     write(unitwf,'(3(1x,e24.17))') (rxyz(j,iat),j=1,3)
!     enddo
!     write(unitwf,*) nvctr_c, nvctr_f
!  else
!     write(unitwf) iorb,eval
!     write(unitwf) hx,hy,hz
!     write(unitwf) n1,n2,n3
!     write(unitwf) locregCenter(1),locregCenter(2),locregCenter(3),locrad,confPotOrder, confPotprefac
!     write(unitwf) nat
!     do iat=1,nat
!     write(unitwf) (rxyz(j,iat),j=1,3)
!     enddo
!     write(unitwf) nvctr_c, nvctr_f
!  end if
!
!  ! coarse part
!  do iseg=1,nseg_c
!     jj=keyv_c(iseg)
!     j0=keyg_c(1,iseg)
!     j1=keyg_c(2,iseg)
!     ii=j0-1
!     i3=ii/((n1+1)*(n2+1))
!     ii=ii-i3*(n1+1)*(n2+1)
!     i2=ii/(n1+1)
!     i0=ii-i2*(n1+1)
!     i1=i0+j1-j0
!     do i=i0,i1
!        tt=psi_c(i-i0+jj)
!        if (useFormattedOutput) then
!           write(unitwf,'(3(i4),1x,e19.12)') i,i2,i3,tt
!        else
!           write(unitwf) i,i2,i3,tt
!        end if
!     enddo
!  enddo
!
!  ! fine part
!  do iseg=1,nseg_f
!     jj=keyv_f(iseg)
!     j0=keyg_f(1,iseg)
!     j1=keyg_f(2,iseg)
!     ii=j0-1
!     i3=ii/((n1+1)*(n2+1))
!     ii=ii-i3*(n1+1)*(n2+1)
!     i2=ii/(n1+1)
!     i0=ii-i2*(n1+1)
!     i1=i0+j1-j0
!     do i=i0,i1
!        t1=psi_f(1,i-i0+jj)
!        t2=psi_f(2,i-i0+jj)
!        t3=psi_f(3,i-i0+jj)
!        t4=psi_f(4,i-i0+jj)
!        t5=psi_f(5,i-i0+jj)
!        t6=psi_f(6,i-i0+jj)
!        t7=psi_f(7,i-i0+jj)
!        if (useFormattedOutput) then
!           write(unitwf,'(3(i4),7(1x,e17.10))') i,i2,i3,t1,t2,t3,t4,t5,t6,t7
!        else
!           write(unitwf) i,i2,i3,t1,t2,t3,t4,t5,t6,t7
!        end if
!     enddo
!  enddo
!
!  if (verbose >= 2) write(*,'(1x,i0,a)') iorb,'th wavefunction written'
!
!END SUBROUTINE writeonewave_linear
!
!subroutine writeLinearCoefficients(unitwf,useFormattedOutput,n1,n2,n3,hx,hy,hz,nat,rxyz,&
!           norb,ntmb,nvctr_c,nvctr_f,coeff)
!  use module_base
!  implicit none
!  logical, intent(in) :: useFormattedOutput
!  integer, intent(in) :: unitwf,norb,n1,n2,n3,nat,ntmb,nvctr_c,nvctr_f
!  real(gp), intent(in) :: hx,hy,hz
!  real(wp), dimension(ntmb,norb), intent(in) :: coeff
!  real(gp), dimension(3,nat), intent(in) :: rxyz
!  !local variables
!  integer :: iat,jj,j0,j1,ii,i0,i1,i2,i3,i,iseg,j
!  real(wp) :: tt,t1,t2,t3,t4,t5,t6,t7
!
!  ! Write the Header
!  if (useFormattedOutput) then
!     write(unitwf,*) norb,ntmb
!     write(unitwf,*) hx,hy,hz
!     write(unitwf,*) n1,n2,n3
!     write(unitwf,*) nat
!     do iat=1,nat
!     write(unitwf,'(3(1x,e24.17))') (rxyz(j,iat),j=1,3)
!     enddo
!     write(unitwf,*) nvctr_c, nvctr_f
!  else
!     write(unitwf) norb, ntmb
!     write(unitwf) hx,hy,hz
!     write(unitwf) n1,n2,n3
!     write(unitwf) nat
!     do iat=1,nat
!     write(unitwf) (rxyz(j,iat),j=1,3)
!     enddo
!     write(unitwf) nvctr_c, nvctr_f
!  end if
!
!  ! Now write the coefficients
!  do i = 1, norb
!     do j = 1, ntmb
!        tt = coeff(j,i)
!        if (useFormattedOutput) then
!           write(unitwf,'(2(i4),1x,e19.12)') i,j,tt
!        else
!           write(unitwf) i,j,tt
!        end if
!     end do
!  end do  
!
!  if (verbose >= 2) write(*,'(1x,a)') 'Wavefunction coefficients written'
!
!END SUBROUTINE writeLinearCoefficients
!
!subroutine writemywaves_linear(iproc,filename,iformat,Lzd,orbs,norb,hx,hy,hz,at,rxyz,psi,coeff)
!  use module_types
!  use module_base
!  use module_interfaces, except_this_one => writeonewave
!  implicit none
!  integer, intent(in) :: iproc,iformat
!  integer, intent(in) :: norb   
!  real(gp), intent(in) :: hx,hy,hz
!  type(atoms_data), intent(in) :: at
!  type(orbitals_data), intent(in) :: orbs         
!  type(local_zone_descriptors), intent(in) :: Lzd
!  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
!  real(wp), dimension(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f,orbs%nspinor,orbs%norbp), intent(in) :: psi  ! Should be the real linear dimension and not the global
!  real(wp), dimension(orbs%norb,norb), intent(in) :: coeff
!  character(len=*), intent(in) :: filename
!  !Local variables
!  integer :: ncount1,ncount_rate,ncount_max,iorb,ncount2,iorb_out,ispinor,ilr
!  real(kind=4) :: tr0,tr1
!  real(kind=8) :: tel
!
!  if (iproc == 0) write(*,"(1x,A,A,a)") "Write wavefunctions to file: ", trim(filename),'.*'
!
!  if (iformat == WF_FORMAT_ETSF) then
!      stop 'Linear scaling with ETSF writing not implemented yet'
!     call write_waves_etsf(iproc,filename,orbs,n1,n2,n3,hx,hy,hz,at,rxyz,wfd,psi)
!  else
!     call cpu_time(tr0)
!     call system_clock(ncount1,ncount_rate,ncount_max)
!
!     ! Write the TMBs in the Plain BigDFT files.
!     do iorb=1,orbs%norbp
!        ilr = orbs%inwhichlocreg(iorb+orbs%isorb)
!        do ispinor=1,orbs%nspinor
!           call open_filename_of_iorb(99,(iformat == WF_FORMAT_BINARY),filename, &
!                & orbs,iorb,ispinor,iorb_out)
!           call writeonewave_linear(99,(iformat == WF_FORMAT_PLAIN),iorb_out,Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3,&
!                Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),Lzd%Llr(ilr)%locregCenter,Lzd%Llr(ilr)%locrad, 4, 0.0d0, &  !put here the real potentialPrefac and Order
!                at%astruct%nat,rxyz,Lzd%Llr(ilr)%wfd%nseg_c,Lzd%Llr(ilr)%wfd%nvctr_c,&
!                Lzd%Llr(ilr)%wfd%keyglob(1,1),Lzd%Llr(ilr)%wfd%keyvglob(1),Lzd%Llr(ilr)%wfd%nseg_f,Lzd%Llr(ilr)%wfd%nvctr_f,&
!                Lzd%Llr(ilr)%wfd%keyglob(1,Lzd%Llr(ilr)%wfd%nseg_c+1),Lzd%Llr(ilr)%wfd%keyvglob(Lzd%Llr(ilr)%wfd%nseg_c+1), &
!                psi(1,ispinor,iorb),psi(Lzd%Llr(ilr)%wfd%nvctr_c+1,ispinor,iorb),orbs%eval(iorb+orbs%isorb))
!           close(unit=99)
!        end do
!     enddo
!
!    ! Now write the coefficients to file
!    ! Must be careful, the orbs%norb is the number of basis functions
!    ! while the norb is the number of orbitals.
!    if(iproc == 0) then
!      call writeLinearCoefficients(99,(iformat == WF_FORMAT_PLAIN),Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3,&
!           Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),at%astruct%nat,rxyz,norb,orbs%norb,Lzd%Glr%wfd%nvctr_c,Lzd%Glr%wfd%nvctr_f,coeff)
!    end if
!     call cpu_time(tr1)
!     call system_clock(ncount2,ncount_rate,ncount_max)
!     tel=dble(ncount2-ncount1)/dble(ncount_rate)
!     write(*,'(a,i4,2(1x,1pe10.3))') '- WRITE WAVES TIME',iproc,tr1-tr0,tel
!     !write(*,'(a,1x,i0,a)') '- iproc',iproc,' finished writing waves'
!  end if
!
!END SUBROUTINE writemywaves_linear
!
!subroutine readonewave_linear(unitwf,useFormattedInput,iorb,iproc,n1,n2,n3,&
!     & hx,hy,hz,at,wfd,rxyz_old,rxyz,locrad,locregCenter,confPotOrder,&
!     & confPotprefac,psi,eval,psifscf)
!  use module_base
!  use module_types
!  use internal_io
!  use module_interfaces
!  implicit none
!  logical, intent(in) :: useFormattedInput
!  integer, intent(in) :: unitwf,iorb,iproc,n1,n2,n3
!  type(wavefunctions_descriptors), intent(in) :: wfd
!  type(atoms_data), intent(in) :: at
!  real(gp), intent(in) :: hx,hy,hz
!  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
!  integer, intent(out) :: confPotOrder
!  real(gp), intent(out) :: locrad, confPotprefac
!  real(wp), intent(out) :: eval
!  real(gp), dimension(3), intent(out) :: locregCenter
!  real(gp), dimension(3,at%astruct%nat), intent(out) :: rxyz_old
!  real(wp), dimension(wfd%nvctr_c+7*wfd%nvctr_f), intent(out) :: psi
!  real(wp), dimension(*), intent(out) :: psifscf !this supports different BC
!  
!  !local variables
!  character(len=*), parameter :: subname='readonewave_linear'
!  character(len = 256) :: error
!  logical :: perx,pery,perz,lstat
!  integer :: iorb_old,n1_old,n2_old,n3_old,iat,iel,nvctr_c_old,nvctr_f_old,i_stat,i_all,i1,i2,i3
!  real(wp) :: tt,t1,t2,t3,t4,t5,t6,t7
!  real(gp) :: tx,ty,tz,displ,hx_old,hy_old,hz_old,mindist
!  real(wp), dimension(:,:,:,:,:,:), allocatable :: psigold
!
!  !write(*,*) 'INSIDE readonewave'
!  call io_read_descr_linear(unitwf, useFormattedInput, iorb_old, eval, n1_old, n2_old, n3_old, &
!       & hx_old, hy_old, hz_old, lstat, error, nvctr_c_old, nvctr_f_old, rxyz_old, at%astruct%nat,&
!       & locrad, locregCenter, confPotOrder, confPotprefac)
!  if (.not. lstat) call io_error(trim(error))
!  if (iorb_old /= iorb) stop 'readonewave_linear'
!
!  !conditions for periodicity in the three directions
!  perx=(at%geocode /= 'F')
!  pery=(at%geocode == 'P')
!  perz=(at%geocode /= 'F')
!
!  tx=0.0_gp
!  ty=0.0_gp
!  tz=0.0_gp
!  do iat=1,at%astruct%nat
!     tx=tx+mindist(perx,at%astruct%cell_dim(1),rxyz(1,iat),rxyz_old(1,iat))**2
!     ty=ty+mindist(pery,at%astruct%cell_dim(2),rxyz(2,iat),rxyz_old(2,iat))**2
!     tz=tz+mindist(perz,at%astruct%cell_dim(3),rxyz(3,iat),rxyz_old(3,iat))**2
!  enddo
!  displ=sqrt(tx+ty+tz)
!
!  if (hx_old == hx .and. hy_old == hy .and. hz_old == hz .and.&
!       n1_old == n1  .and. n2_old == n2 .and. n3_old == n3 .and. displ <= 1.d-3) then
!
!     if (iproc == 0) write(*,*) 'wavefunctions need NO reformatting'
!     call read_psi_compress(unitwf, useFormattedInput, nvctr_c_old, nvctr_f_old, psi, lstat, error)
!     if (.not. lstat) call io_error(trim(error))
!
!  else
!
!     if (iproc == 0 .and. iorb == 1) then
!        write(*,*) 'wavefunctions need reformatting'
!        if (hx_old /= hx .or. hy_old /= hy .or. hz_old /= hz) write(*,"(1x,A,6F14.10)") &
!             'because hgrid_old /= hgrid',hx_old,hy_old,hz_old,hx,hy,hz
!        if (n1_old /= n1  .or. n2_old /= n2 .or. n3_old /= n3 ) &
!             write(*,*) 'because cell size has changed',n1_old,n1,n2_old,n2,n3_old,n3
!        if (displ > 1.d-3 ) write(*,*) 'large displacement of molecule',displ
!     end if
!
! NOT SURE YET WHAT SHOULD BE DONE FOR LINEAR CASE, so just stop
!if (iproc == 0) write(*,*) 'This is forbiden for now in linear case!'
!call mpi_finalize(i_all)
!stop 
!
!     allocate(psigold(0:n1_old,2,0:n2_old,2,0:n3_old,2+ndebug),stat=i_stat)
!     call memocc(i_stat,psigold,'psigold',subname)
!
!     call to_zero(8*(n1_old+1)*(n2_old+1)*(n3_old+1),psigold)
!     do iel=1,nvctr_c_old
!        if (useFormattedInput) then
!           read(unitwf,*) i1,i2,i3,tt
!        else
!           read(unitwf) i1,i2,i3,tt
!        end if
!        psigold(i1,1,i2,1,i3,1)=tt
!     enddo
!     do iel=1,nvctr_f_old
!        if (useFormattedInput) then
!           read(unitwf,*) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
!        else
!           read(unitwf) i1,i2,i3,t1,t2,t3,t4,t5,t6,t7
!        end if
!        psigold(i1,2,i2,1,i3,1)=t1
!        psigold(i1,1,i2,2,i3,1)=t2
!        psigold(i1,2,i2,2,i3,1)=t3
!        psigold(i1,1,i2,1,i3,2)=t4
!        psigold(i1,2,i2,1,i3,2)=t5
!        psigold(i1,1,i2,2,i3,2)=t6
!        psigold(i1,2,i2,2,i3,2)=t7
!     enddo
!
!     ! I put nat = 1 here, since only one position is saved in wavefunction files.
!     call reformatonewave(displ,wfd,at,hx_old,hy_old,hz_old,n1_old,n2_old,n3_old,&
!          rxyz_old,psigold,hx,hy,hz,n1,n2,n3,rxyz,psifscf,psi)
!
!     i_all=-product(shape(psigold))*kind(psigold)
!     deallocate(psigold,stat=i_stat)
!     call memocc(i_stat,i_all,'psigold',subname)
!
!  endif
!
!END SUBROUTINE readonewave_linear                                                     
!
!subroutine io_read_descr_linear(unitwf, formatted, iorb_old, eval, n1_old, n2_old, n3_old, &
!       & hx_old, hy_old, hz_old, lstat, error, nvctr_c_old, nvctr_f_old, rxyz_old, nat, &
!       & locrad, locregCenter, confPotOrder, confPotprefac)
!    use module_base
!    use module_types
!    use internal_io
!    implicit none
!
!    integer, intent(in) :: unitwf
!    logical, intent(in) :: formatted
!    integer, intent(out) :: iorb_old
!    integer, intent(out) :: n1_old, n2_old, n3_old
!    real(gp), intent(out) :: hx_old, hy_old, hz_old
!    logical, intent(out) :: lstat
!    real(wp), intent(out) :: eval
!    integer, intent(out) :: confPotOrder
!    real(gp), intent(out) :: locrad, confPotprefac
!    real(gp), dimension(3), intent(out) :: locregCenter
!    character(len =256), intent(out) :: error
!    ! Optional arguments
!    integer, intent(out), optional :: nvctr_c_old, nvctr_f_old
!    integer, intent(in), optional :: nat
!    real(gp), dimension(:,:), intent(out), optional :: rxyz_old
!
!    character(len = *), parameter :: subname = "io_read_descr_linear"
!    integer :: i, iat, i_stat, nat_
!    real(gp) :: rxyz(3)
!
!    lstat = .false.
!    write(error, "(A)") "cannot read psi description."
!    if (formatted) then
!       read(unitwf,*,iostat=i_stat) iorb_old,eval
!       if (i_stat /= 0) return
!       read(unitwf,*,iostat=i_stat) hx_old,hy_old,hz_old
!       if (i_stat /= 0) return
!       read(unitwf,*,iostat=i_stat) n1_old,n2_old,n3_old
!       if (i_stat /= 0) return
!       read(unitwf,*,iostat=i_stat) (locregCenter(i),i=1,3),locrad,confPotOrder, confPotprefac
!       if (i_stat /= 0) return
!       !write(*,*) 'reading ',nat,' atomic positions'
!       if (present(nat) .And. present(rxyz_old)) then
!          read(unitwf,*,iostat=i_stat) nat_
!          if (i_stat /= 0) return
!          ! Sanity check
!          if (size(rxyz_old, 2) /= nat) stop "Mismatch in coordinate array size."
!          if (nat_ /= nat) stop "Mismatch in coordinate array size."
!          do iat=1,nat
!             read(unitwf,*,iostat=i_stat) (rxyz_old(i,iat),i=1,3)
!             if (i_stat /= 0) return
!          enddo
!       else
!          read(unitwf,*,iostat=i_stat) nat_
!          if (i_stat /= 0) return
!          do iat=1,nat_
!             read(unitwf,*,iostat=i_stat)
!             if (i_stat /= 0) return
!          enddo
!       end if
!       if (present(nvctr_c_old) .and. present(nvctr_f_old)) then
!          read(unitwf,*,iostat=i_stat) nvctr_c_old, nvctr_f_old
!          if (i_stat /= 0) return
!       else
!          read(unitwf,*,iostat=i_stat) i, iat
!          if (i_stat /= 0) return
!       end if
!    else
!       read(unitwf,iostat=i_stat) iorb_old,eval
!       if (i_stat /= 0) return
!       read(unitwf,iostat=i_stat) hx_old,hy_old,hz_old
!       if (i_stat /= 0) return
!       read(unitwf,iostat=i_stat) n1_old,n2_old,n3_old
!       if (i_stat /= 0) return
!       read(unitwf,iostat=i_stat) (locregCenter(i),i=1,3),locrad,confPotOrder, confPotprefac
!       if (i_stat /= 0) return
!       if (present(nat) .And. present(rxyz_old)) then
!          read(unitwf,iostat=i_stat) nat_
!          if (i_stat /= 0) return
!          ! Sanity check
!          if (size(rxyz_old, 2) /= nat) stop "Mismatch in coordinate array size." 
!          if (nat_ /= nat) stop "Mismatch in coordinate array size."
!          do iat=1,nat
!             read(unitwf,iostat=i_stat)(rxyz_old(i,iat),i=1,3)
!             if (i_stat /= 0) return
!          enddo
!       else
!          read(unitwf,iostat=i_stat) nat_
!          if (i_stat /= 0) return
!          do iat=1,nat_
!             read(unitwf,iostat=i_stat) rxyz
!             if (i_stat /= 0) return
!          enddo
!       end if
!       if (present(nvctr_c_old) .and. present(nvctr_f_old)) then
!          read(unitwf,iostat=i_stat) nvctr_c_old, nvctr_f_old
!          if (i_stat /= 0) return
!       else
!          read(unitwf,iostat=i_stat) i, iat
!          if (i_stat /= 0) return
!       end if
!    end if
!    lstat = .true.
!END SUBROUTINE io_read_descr_linear
!
!subroutine io_read_descr_coeff(unitwf, formatted, norb_old, ntmb_old, n1_old, n2_old, n3_old, &
!       & hx_old, hy_old, hz_old, lstat, error, nvctr_c_old, nvctr_f_old, rxyz_old, nat)
!    use module_base
!    use module_types
!    use internal_io
!    implicit none
!    integer, intent(in) :: unitwf
!    logical, intent(in) :: formatted
!    integer, intent(out) :: norb_old, ntmb_old
!    integer, intent(out) :: n1_old, n2_old, n3_old
!    real(gp), intent(out) :: hx_old, hy_old, hz_old
!    logical, intent(out) :: lstat
!    character(len =256), intent(out) :: error
!    ! Optional arguments
!    integer, intent(out), optional :: nvctr_c_old, nvctr_f_old
!    integer, intent(in), optional :: nat
!    real(gp), dimension(:,:), intent(out), optional :: rxyz_old
!
!    character(len = *), parameter :: subname = "io_read_descr_linear"
!    integer :: i, iat, i_stat, nat_
!    real(gp) :: rxyz(3)
!
!    lstat = .false.
!    write(error, "(A)") "cannot read psi description."
!    if (formatted) then
!       read(unitwf,*,iostat=i_stat) norb_old , ntmb_old
!       if (i_stat /= 0) return
!       read(unitwf,*,iostat=i_stat) hx_old,hy_old,hz_old
!       if (i_stat /= 0) return
!       read(unitwf,*,iostat=i_stat) n1_old,n2_old,n3_old
!       if (i_stat /= 0) return
!       !write(*,*) 'reading ',nat,' atomic positions'
!       if (present(nat) .And. present(rxyz_old)) then
!          read(unitwf,*,iostat=i_stat) nat_
!          if (i_stat /= 0) return
!          ! Sanity check
!          if (size(rxyz_old, 2) /= nat) stop "Mismatch in coordinate array size."
!          if (nat_ /= nat) stop "Mismatch in coordinate array size."
!          do iat=1,nat
!             read(unitwf,*,iostat=i_stat) (rxyz_old(i,iat),i=1,3)
!             if (i_stat /= 0) return
!          enddo
!       else
!          read(unitwf,*,iostat=i_stat) nat_
!          if (i_stat /= 0) return
!          do iat=1,nat_
!             read(unitwf,*,iostat=i_stat)
!             if (i_stat /= 0) return
!          enddo
!       end if
!       if (present(nvctr_c_old) .and. present(nvctr_f_old)) then
!          read(unitwf,*,iostat=i_stat) nvctr_c_old, nvctr_f_old
!          if (i_stat /= 0) return
!       else
!          read(unitwf,*,iostat=i_stat) i, iat
!          if (i_stat /= 0) return
!       end if
!    else
!       read(unitwf,iostat=i_stat) norb_old, ntmb_old
!       if (i_stat /= 0) return
!       read(unitwf,iostat=i_stat) hx_old,hy_old,hz_old
!       if (i_stat /= 0) return
!       read(unitwf,iostat=i_stat) n1_old,n2_old,n3_old
!       if (i_stat /= 0) return
!       if (present(nat) .And. present(rxyz_old)) then
!          read(unitwf,iostat=i_stat) nat_
!          if (i_stat /= 0) return
!          ! Sanity check
!          if (size(rxyz_old, 2) /= nat) stop "Mismatch in coordinate array size." 
!          if (nat_ /= nat) stop "Mismatch in coordinate array size."
!          do iat=1,nat
!             read(unitwf,iostat=i_stat)(rxyz_old(i,iat),i=1,3)
!             if (i_stat /= 0) return
!          enddo
!       else
!          read(unitwf,iostat=i_stat) nat_
!          if (i_stat /= 0) return
!          do iat=1,nat_
!             read(unitwf,iostat=i_stat) rxyz
!             if (i_stat /= 0) return
!          enddo
!       end if
!       if (present(nvctr_c_old) .and. present(nvctr_f_old)) then
!          read(unitwf,iostat=i_stat) nvctr_c_old, nvctr_f_old
!          if (i_stat /= 0) return
!       else
!          read(unitwf,iostat=i_stat) i, iat
!          if (i_stat /= 0) return
!       end if
!    end if
!    lstat = .true.
!END SUBROUTINE io_read_descr_coeff
!
!
!subroutine read_coeff_minbasis(unitwf,useFormattedInput,iproc,n1,n2,n3,norb,ntmb,&
!     & hx,hy,hz,at,rxyz_old,rxyz,coeff)
!  use module_base
!  use module_types
!  use internal_io
!  use module_interfaces
!  implicit none
!  logical, intent(in) :: useFormattedInput
!  integer, intent(in) :: unitwf,iproc,n1,n2,n3,norb,ntmb
!  type(atoms_data), intent(in) :: at
!  real(gp), intent(in) :: hx,hy,hz
!  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
!  real(gp), dimension(3,at%astruct%nat), intent(out) :: rxyz_old
!  real(wp), dimension(ntmb,norb), intent(out) :: coeff
!  !local variables
!  character(len=*), parameter :: subname='readonewave_linear'
!  character(len = 256) :: error
!  logical :: perx,pery,perz,lstat
!  integer :: norb_old,n1_old,n2_old,n3_old,iat,nvctr_c_old,nvctr_f_old,i_stat,i_all
!  integer :: ntmb_old, i1, i2,i,j
!  real(wp) :: tt
!  real(gp) :: tx,ty,tz,displ,hx_old,hy_old,hz_old,mindist
!
!  !write(*,*) 'INSIDE readonewave'
!  call io_read_descr_coeff(unitwf, useFormattedInput, norb_old, ntmb_old, n1_old, n2_old, n3_old, &
!       & hx_old, hy_old, hz_old, lstat, error, nvctr_c_old, nvctr_f_old, rxyz_old, at%astruct%nat)
!  if (.not. lstat) call io_error(trim(error))
!
!  !conditions for periodicity in the three directions
!  perx=(at%geocode /= 'F')
!  pery=(at%geocode == 'P')
!  perz=(at%geocode /= 'F')
!
!  tx=0.0_gp
!  ty=0.0_gp
!  tz=0.0_gp
!  do iat=1,at%astruct%nat
!     tx=tx+mindist(perx,at%astruct%cell_dim(1),rxyz(1,iat),rxyz_old(1,iat))**2
!     ty=ty+mindist(pery,at%astruct%cell_dim(2),rxyz(2,iat),rxyz_old(2,iat))**2
!     tz=tz+mindist(perz,at%astruct%cell_dim(3),rxyz(3,iat),rxyz_old(3,iat))**2
!  enddo
!  displ=sqrt(tx+ty+tz)
!
!  if (hx_old == hx .and. hy_old == hy .and. hz_old == hz .and.&
!       n1_old == n1  .and. n2_old == n2 .and. n3_old == n3 .and. displ <= 1.d-3 .and. &
!       norb == norb_old .and. ntmb == ntmb_old) then
!
!     if (iproc == 0) write(*,*) 'wavefunctions need NO reformatting'
!
!     ! Now write the coefficients
!     do i = 1, norb
!        do j = 1, ntmb
!           if (useFormattedInput) then
!              read(unitwf,*,iostat=i_stat) i1,i2,tt
!           else
!              read(unitwf,iostat=i_stat) i1,i2,tt
!           end if
!           if (i_stat /= 0) stop 'Problem reading the coefficients'
!           coeff(j,i) = tt  
!        end do
!     end do
!     if (verbose >= 2) write(*,'(1x,a)') 'Wavefunction coefficients written'
!
!  else
!     if (iproc == 0) then
!        write(*,*) 'wavefunctions need reformatting'
!        if (hx_old /= hx .or. hy_old /= hy .or. hz_old /= hz) write(*,"(1x,A,6F14.10)") &
!             'because hgrid_old /= hgrid',hx_old,hy_old,hz_old,hx,hy,hz
!        if (n1_old /= n1  .or. n2_old /= n2 .or. n3_old /= n3 ) &
!             write(*,*) 'because cell size has changed',n1_old,n1,n2_old,n2,n3_old,n3
!        if (displ > 1.d-3 ) write(*,*) 'large displacement of molecule',displ
!        if (norb /= norb_old) write(*,*) 'Differing number of orbitals',norb,norb_old
!        if (ntmb /= ntmb_old) write(*,*) 'Differing number of minimal basis functions',ntmb,ntmb_old
!     end if
!
!     ! NOT SURE YET WHAT SHOULD BE DONE FOR LINEAR CASE, so just stop
!     if (iproc == 0) then
!        write(*,*) 'This is forbiden for now in linear case!'
!        call mpi_finalize(i_all)
!        stop
!     end if
!  end if
!
!END SUBROUTINE read_coeff_minbasis
!
!
!  Reads wavefunction from file and transforms it properly if hgrid or size of simulation cell
!  have changed
!subroutine readmywaves_linear(iproc,filename,iformat,norb,Lzd,orbs,at,rxyz_old,rxyz,  & 
!    psi,coeff,orblist)
!  use module_base
!  use module_types
!  use module_interfaces, except_this_one => readmywaves_linear
!  implicit none
!  integer, intent(in) :: iproc, iformat,norb
!  type(orbitals_data), intent(inout) :: orbs  ! orbs related to the basis functions
!  type(local_zone_descriptors), intent(in) :: Lzd
!  type(atoms_data), intent(in) :: at
!  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
!  real(gp), dimension(3,at%astruct%nat), intent(out) :: rxyz_old
!  real(wp), dimension(orbs%npsidim_orbs), intent(out) :: psi  
!  real(gp), dimension(norb,orbs%norb),intent(out) :: coeff
!  character(len=*), intent(in) :: filename
!  integer, dimension(orbs%norb), optional :: orblist
!  !Local variables
!  character(len=*), parameter :: subname='readmywaves_linear'
!  logical :: perx,pery,perz
!  integer :: ncount1,ncount_rate,ncount_max,iorb,i_stat,i_all,ncount2,nb1,nb2,nb3
!  integer :: iorb_out,ispinor,ilr,ind
!  integer :: confPotOrder
!  real(gp) :: locrad, confPotprefac
!  real(gp), dimension(3) :: locregCenter
!  real(kind=4) :: tr0,tr1
!  real(kind=8) :: tel
!  real(wp), dimension(:,:,:), allocatable :: psifscf
!  !integer, dimension(orbs%norb) :: orblist2
!
!  call cpu_time(tr0)
!  call system_clock(ncount1,ncount_rate,ncount_max)
!
!  if (iformat == WF_FORMAT_ETSF) then
!     stop 'Linear scaling with ETSF writing not implemented yet'
!     !construct the orblist or use the one in argument
!     !do nb1 = 1, orbs%norb
!     !orblist2(nb1) = nb1
!     !if(present(orblist)) orblist2(nb1) = orblist(nb1) 
!     !end do
!
!     !call read_waves_etsf(iproc,filename // ".etsf",orbs,n1,n2,n3,hx,hy,hz,at,rxyz_old,rxyz,  & 
!     !     wfd,psi)
!  else if (iformat == WF_FORMAT_BINARY .or. iformat == WF_FORMAT_PLAIN) then
!     !conditions for periodicity in the three directions
!     !perx=(at%geocode /= 'F')
!     !pery=(at%geocode == 'P')
!     !perz=(at%geocode /= 'F')
!
!     !buffers related to periodicity
!     !WARNING: the boundary conditions are not assumed to change between new and old
!     !call ext_buffers_coarse(perx,nb1)
!     !call ext_buffers_coarse(pery,nb2)
!     !call ext_buffers_coarse(perz,nb3)
!     !allocate(psifscf(-nb1:2*n1+1+nb1,-nb2:2*n2+1+nb2,-nb3:2*n3+1+nb3+ndebug),stat=i_stat)
!     !call memocc(i_stat,psifscf,'psifscf',subname)
!     allocate(psifscf(1,1,1+ndebug),stat=i_stat)
!     call memocc(i_stat,psifscf,'psifscf',subname)
!     ind = 0
!     do iorb=1,orbs%norbp!*orbs%nspinor
!        ilr = orbs%inwhichlocreg(iorb+orbs%isorb)
!        do ispinor=1,orbs%nspinor
!           if(present(orblist)) then
!              call open_filename_of_iorb(99,(iformat == WF_FORMAT_BINARY),filename, &
!                   & orbs,iorb,ispinor,iorb_out, orblist(iorb+orbs%isorb))
!           else
!              call open_filename_of_iorb(99,(iformat == WF_FORMAT_BINARY),filename, &
!                   & orbs,iorb,ispinor,iorb_out)
!           end if           
!           call readonewave_linear(99, (iformat == WF_FORMAT_PLAIN),iorb_out,iproc,&
!                Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3,Lzd%hgrids(1),Lzd%hgrids(2),&
!                Lzd%hgrids(3),at,Lzd%Llr(ilr)%wfd,rxyz_old,rxyz,locrad,locregCenter,&
!                confPotOrder,confPotPrefac,psi(1+ind),orbs%eval(orbs%isorb+iorb),psifscf)
!           close(unit=99)
!           ind = ind + Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f
!        end do
!
!     end do
!
!     i_all=-product(shape(psifscf))*kind(psifscf)
!     deallocate(psifscf,stat=i_stat)
!     call memocc(i_stat,i_all,'psifscf',subname)
!
!     !Open the coefficient file 
!     if(iformat == WF_FORMAT_PLAIN) then
!        open(unit=99,file=filename//'_coeff.bin',status='unknown',form='formatted')
!     else if(iformat == WF_FORMAT_BINARY) then
!        open(unit=99,file=filename//'_coeff.bin',status='unknown',form='unformatted')
!     else
!        stop 'Coefficient format not implemented'
!     end if
!write(*,*) 'arriving here'
!     call read_coeff_minbasis(99,(iformat == WF_FORMAT_PLAIN),iproc,Lzd%Glr%d%n1,Lzd%Glr%d%n2,Lzd%Glr%d%n3,norb,orbs%norb,&
!     & Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3),at,rxyz_old,rxyz,coeff)
!     close(unit=99)
!  else
!     write(0,*) "Unknown wavefunction file format from filename."
!     stop
!  end if
!
!  call cpu_time(tr1)
!  call system_clock(ncount2,ncount_rate,ncount_max)
!  tel=dble(ncount2-ncount1)/dble(ncount_rate)
!  write(*,'(a,i4,2(1x,1pe10.3))') '- READING WAVES TIME',iproc,tr1-tr0,tel
!END SUBROUTINE readmywaves_linear
!
!
!subroutine initialize_linear_from_file(iproc,nproc,filename,iformat,Lzd,orbs,at,rxyz,orblist)
!  use module_base
!  use module_types
!  use module_defs
!  use module_interfaces, except_this_one => initialize_linear_from_file
!  implicit none
!  integer, intent(in) :: iproc, nproc, iformat
!  type(orbitals_data), intent(inout) :: orbs  
!  type(atoms_data), intent(in) :: at
!  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
!  character(len=*), intent(in) :: filename
!  type(local_zone_descriptors), intent(inout) :: Lzd 
!  integer, dimension(orbs%norb), optional :: orblist
!  !Local variables
!  character(len=*), parameter :: subname='initialize_linear_from_file'
!  character(len =256) :: error
!  logical :: lstat, consistent, perx, pery, perz
!  integer :: ilr, ierr, iorb_old, iorb, jorb, ispinor, iorb_out, n1_old, n2_old, n3_old
!  integer :: nlr, nvctr_c_old, nvctr_f_old, i_stat, i_all,confPotOrder, confPotOrder_old
!  real(kind=4) :: tr0,tr1
!  real(kind=8) :: tel,dx,dy,dz,dist,eval
!  real(gp) :: hx_old, hy_old, hz_old, mindist
!  real(gp), dimension(orbs%norb):: locrad, confPotprefac
!  real(gp), dimension(3,at%astruct%nat) :: rxyz_old
!  real(gp), dimension(3,orbs%norb) :: locregCenter
!  integer, dimension(:), allocatable :: lrtable
!  integer, dimension(orbs%norb) :: nvctr_c, nvctr_f
!  real(gp), dimension(:), allocatable :: lrad
!  real(gp), dimension(:,:), allocatable :: cxyz
!  logical, dimension(:), allocatable :: calcbounds
!
!  ! NOTES:
!  ! The orbs%norb family must be all constructed before this routine
!  ! This can be done from the input.lin since the number of basis functions should be fixed.
!
!  call to_zero(3*orbs%norb,locregCenter(1,1))
!  call to_zero(orbs%norb,locrad(1))
!  call to_zero(orbs%norb,confPotprefac(1))
!
!  ! First read the headers (reading is distributed) and then the information is communicated to all procs.
!  ! Then each proc generates a group of lrs that are communicated to all others.
!  if (iformat == WF_FORMAT_ETSF) then
!     stop 'Linear scaling with ETSF writing not implemented yet'
!  else if (iformat == WF_FORMAT_BINARY .or. iformat == WF_FORMAT_PLAIN) then
!     do iorb=1,orbs%norbp!*orbs%nspinor
!        do ispinor=1,orbs%nspinor
!           if(present(orblist)) then
!              call open_filename_of_iorb(99,(iformat == WF_FORMAT_BINARY),filename, &
!                   & orbs,iorb,ispinor,iorb_out, orblist(iorb+orbs%isorb))
!           else
!              call open_filename_of_iorb(99,(iformat == WF_FORMAT_BINARY),filename, &
!                   & orbs,iorb,ispinor,iorb_out)
!           end if          
!           call io_read_descr_linear(99,(iformat == WF_FORMAT_PLAIN), iorb_old, eval, n1_old, n2_old, n3_old, &
!                & hx_old, hy_old, hz_old, lstat, error, nvctr_c(iorb+orbs%isorb), nvctr_f(iorb+orbs%isorb),&
!                & rxyz_old, at%astruct%nat, locrad(iorb+orbs%isorb), locregCenter(1,iorb+orbs%isorb), confPotOrder,&
!                & confPotprefac(iorb+orbs%isorb))
!           if (.not. lstat) then ; write(*,*) trim(error) ; stop; end if
!           if (iorb_old /= iorb_out) stop 'initialize_linear_from_file'
!           close(unit=99)
!TO DO: confPotOrder_old should be read from input.lin
!           if(iorb==1) confPotOrder_old = confPotOrder
!           call check_consistency(Lzd, at, hx_old, hy_old, hz_old, n1_old, n2_old, n3_old, &
!                rxyz_old,rxyz,confPotOrder,confPotOrder_old,consistent)
!           if(.not. consistent) then
!             write(*,*) 'Inconsistency in file, iorb=',iorb_out
!             call mpi_finalize(ierr)
!             stop
!           end if
!           confPotOrder_old = confPotOrder
!        end do
!     end do
!  else
!     write(0,*) "Unknown wavefunction file format from filename."
!     stop
!  end if
!
!  ! Communication of the quantities
!   call mpiallred(locregCenter(1,1),3*orbs%norb,MPI_SUM,MPI_COMM_WORLD,ierr)
!   call mpiallred(locrad(1),orbs%norb,MPI_SUM,MPI_COMM_WORLD,ierr)
!   call mpiallred(confPotprefac(1),orbs%norb,MPI_SUM,MPI_COMM_WORLD,ierr)
!
!  ! Now that each processor has all the information, we can build the locregs
!  ! Find the number of inequivalent locregs
!  allocate(lrtable(orbs%norb),stat=i_stat)
!  call memocc(i_stat,ilr,'ilr',subname)
!  allocate(orbs%inwhichlocreg(orbs%norb),stat=i_stat)
!  call memocc(i_stat,orbs%inwhichlocreg,'orbs%inwhichlocreg',subname)
!
!  nlr = 0
!  lrtable = 0
!  outer_loop: do iorb = 1, orbs%norb
!     do jorb = iorb+1, orbs%norb
!        dx=mindist(perx,at%astruct%cell_dim(1),locregCenter(1,iorb),locregCenter(1,jorb))**2
!        dy=mindist(pery,at%astruct%cell_dim(2),locregCenter(2,iorb),locregCenter(2,jorb))**2
!        dz=mindist(perz,at%astruct%cell_dim(3),locregCenter(3,iorb),locregCenter(3,jorb))**2
!        dist=sqrt(dx+dy+dz)
!        if(dist < 1.0d-3 .and. abs(locrad(iorb)-locrad(jorb)) < 1.0d-3 .and. &
!           confPotprefac(iorb) == confPotprefac(jorb)) then
!           cycle outer_loop
!        end if
!     end do
!     nlr = nlr + 1
!     lrtable(nlr) = iorb
!  end do outer_loop
!
!  Lzd%nlr = nlr
!  allocate(Lzd%Llr(nlr),stat=i_stat)
!  allocate(lrad(nlr),stat=i_stat)
!  call memocc(i_stat,lrad,'lrad',subname)
!  allocate(cxyz(3,nlr),stat=i_stat)
!  call memocc(i_stat,cxyz,'cxyz',subname)
!  allocate(calcbounds(nlr),stat=i_stat)
!  call memocc(i_stat,calcbounds,'calcbounds',subname)
!  
!  
!  do ilr=1,nlr
!     iorb = lrtable(ilr)
!     lrad(ilr) = locrad(iorb)
!     cxyz(1,ilr) = locregCenter(1,iorb)
!     cxyz(2,ilr) = locregCenter(2,iorb)
!     cxyz(3,ilr) = locregCenter(3,iorb)
!     calcbounds(ilr) = .true.
!     do jorb = 1, orbs%norb
!        dx=mindist(perx,at%astruct%cell_dim(1),locregCenter(1,iorb),locregCenter(1,jorb))**2
!        dy=mindist(pery,at%astruct%cell_dim(2),locregCenter(2,iorb),locregCenter(2,jorb))**2
!        dz=mindist(perz,at%astruct%cell_dim(3),locregCenter(3,iorb),locregCenter(3,jorb))**2
!        dist=sqrt(dx+dy+dz)
!        if(dist < 1.0d-3 .and. abs(locrad(iorb)-locrad(jorb)) < 1.0d-3 .and. &
!           confPotprefac(iorb) == confPotprefac(jorb)) then
!           orbs%inwhichlocreg(jorb) = ilr
!        end if
!     end do
!  end do
!
!  i_all = -product(shape(lrtable))*kind(lrtable)
!  deallocate(lrtable,stat=i_stat)
!  call memocc(i_stat,i_all,'lrtable',subname)
!
!TO DO: CUBIC LOCREGS
!  call determine_locregSphere_parallel(iproc,nproc,Lzd%nlr,cxyz,lrad,Lzd%hgrids(1),&
!       Lzd%hgrids(2),Lzd%hgrids(3),Lzd%Glr,Lzd%Llr,calcbounds) 
!  
!END SUBROUTINE initialize_linear_from_file
!
!subroutine check_consistency(Lzd, at, hx_old, hy_old, hz_old, n1_old, n2_old, n3_old, &
!           rxyz_old,rxyz,confPotOrder,confPotOrder_old,consistent)
!  use module_base
!  use module_types
!  implicit none
!  integer, intent(in) :: confPotOrder,confPotOrder_old, n1_old, n2_old, n3_old
!  type(atoms_data), intent(in) :: at
!  real(gp), intent(in) :: hx_old, hy_old, hz_old
!  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz, rxyz_old
!  type(local_zone_descriptors), intent(in) :: Lzd 
!  logical, intent(out) :: consistent
!  ! Local variables
!  logical :: perx, pery, perz
!  integer :: iat
!  real(gp):: tx, ty, tz, displ, mindist  
!
!  !conditions for periodicity in the three directions
!  perx=(at%geocode /= 'F')
!  pery=(at%geocode == 'P')
!  perz=(at%geocode /= 'F')
!
!  tx=0.0_gp
!  ty=0.0_gp
!  tz=0.0_gp
!  do iat=1,at%astruct%nat
!     tx=tx+mindist(perx,at%astruct%cell_dim(1),rxyz(1,iat),rxyz_old(1,iat))**2
!     ty=ty+mindist(pery,at%astruct%cell_dim(2),rxyz(2,iat),rxyz_old(2,iat))**2
!     tz=tz+mindist(perz,at%astruct%cell_dim(3),rxyz(3,iat),rxyz_old(3,iat))**2
!  enddo
!  displ=sqrt(tx+ty+tz)
!  consistent = .true.
!  if(hx_old /= Lzd%hgrids(1) .or. hy_old /= Lzd%hgrids(2) .or. hz_old /= Lzd%hgrids(3)) then
!    write(*,"(1x,A,6F14.10)") 'Stopping because hgrid_old /= hgrid',hx_old,hy_old,hz_old,&
!         Lzd%hgrids(1),Lzd%hgrids(2),Lzd%hgrids(3)
!    consistent = .false.
!  else if (n1_old /= Lzd%Glr%d%n1  .or. n2_old /= Lzd%Glr%d%n2 .or. n3_old /= Lzd%Glr%d%n3 ) then
!    write(*,"(1x,A,6F14.10)") 'Stopping because global cell size',&
!    n1_old,Lzd%Glr%d%n1,n2_old,Lzd%Glr%d%n2,n3_old,Lzd%Glr%d%n3
!    consistent = .false.
!  else if(displ > 1.d-3 ) then
!    write(*,*) 'Stopping because of large displacement of molecule',displ
!    consistent = .false.
!  else if(confpotOrder /= confPotOrder_old) then
!    write(*,*) 'Stopping because of inconsistent confPotOrder',confPotOrder,confPotOrder_old 
!    consistent = .false.
!  end if
!
!END SUBROUTINE check_consistency
