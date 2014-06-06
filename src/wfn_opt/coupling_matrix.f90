!> @file
!!  Routines related to coupling matrix (TD-DFT Casida's formalism)
!! @author
!!    Copyright (C) 2009-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

subroutine center_of_charge(at,rxyz,cc)
  use module_base
  use module_types
  implicit none
  type(atoms_data), intent(in) :: at
  real(gp), dimension(3,at%astruct%nat), intent(in) :: rxyz
  real(gp), dimension(3), intent(out) :: cc
  !local variables
  integer :: iat,ityp
  real(gp) :: zatom,rx,ry,rz,qtot

  cc(1)=0.0_gp
  cc(2)=0.0_gp
  cc(3)=0.0_gp
  qtot=0.0_gp
  do iat=1,at%astruct%nat
     ityp=at%astruct%iatype(iat)
     zatom=real(at%nelpsp(ityp),gp)
     qtot=qtot+zatom
     !coordinates of the center
     rx=rxyz(1,iat)
     ry=rxyz(2,iat)
     rz=rxyz(3,iat)
     cc(1)=cc(1)+rx*zatom
     cc(2)=cc(2)+ry*zatom
     cc(3)=cc(3)+rz*zatom
  end do
  if (qtot /= 0.0_gp) then
     cc(1)=cc(1)/qtot
     cc(2)=cc(2)/qtot
     cc(3)=cc(3)/qtot
  end if
END SUBROUTINE center_of_charge


!> Calculate the coupling matrix needed for Casida's TDDFT approach
subroutine coupling_matrix_prelim(iproc,nproc,geocode,nspin,lr,orbsocc,orbsvirt,i3s,n3p,&
     hxh,hyh,hzh,chargec,pkernel,dvxcdrho,psirocc,psivirtr)
  use module_base
  use module_types
  use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
  use yaml_output
  implicit none
  character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
  integer, intent(in) :: iproc,nproc,n3p,nspin,i3s
  real(gp), intent(in) :: hxh,hyh,hzh
  real(gp), dimension(3) :: chargec
  type(locreg_descriptors), intent(in) :: lr
  type(orbitals_data), intent(in) :: orbsocc,orbsvirt
  real(wp), dimension(lr%d%n1i*lr%d%n2i*n3p*orbsocc%norb), intent(in) :: psirocc
  real(wp), dimension(lr%d%n1i*lr%d%n2i*n3p*orbsvirt%norb), intent(in) :: psivirtr
  real(wp), dimension(lr%d%n1i,lr%d%n2i,n3p,max((nspin*(nspin+1))/2,2)), intent(in) :: dvxcdrho
  type(coulomb_operator) :: pkernel
  !local variables
  character(len=*), parameter :: subname='coupling_matrix_prelim'
  logical :: tda=.true.,onlyfxc=.false.,dofxc=.true.,perx,pery,perz
  integer :: i_all,i_stat,ierr,imulti,jmulti,jorba,jorbi,index
  integer :: i1,i2,i3p,iorbi,iorba,indi,inda,ind2,ind3,ntda,ispin,jspin
  integer :: ik,jk,nmulti,lwork,info,nbl1,nbl2,nbl3,nbr3,nbr2,nbr1,ndipoles
  real(gp) :: ehart,hfac,x,y,z
  real(wp), dimension(:), allocatable :: omega,work
  real(wp), dimension(:,:), allocatable :: K,Kbig,Kaux,dipoles,fi
  real(wp), dimension(:,:,:), allocatable :: v_ias
  real(wp), dimension(:,:,:,:), allocatable :: rho_ias

  if(iproc==0) call yaml_comment('Linear-Response TDDFT calculations',hfill='-')
  !if(iproc==0) write(*,'(1x,a)')"=========================================================="
  !if(iproc==0) write(*,'(1x,a)')" Linear-Response TDDFT calculations"



  !conditions for periodicity in the three directions
  perx=(geocode /= 'F')
  pery=(geocode == 'P')
  perz=(geocode /= 'F')

  call ext_buffers(perx,nbl1,nbr1)
  call ext_buffers(pery,nbl2,nbr2)
  call ext_buffers(perz,nbl3,nbr3)


  !calculate the number of multiindices
  nmulti=0
  do imulti=1,orbsvirt%norb*orbsocc%norb
     !calculate the orbital index
     iorbi=(imulti-1)/orbsvirt%norb+1
     iorba=imulti-(iorbi-1)*orbsvirt%norb
     !check the spin of the orbitals considered
     if (orbsocc%spinsgn(iorbi) == orbsvirt%spinsgn(iorba)) then
        nmulti=nmulti+1
     end if
  end do

  !partial densities and potentials
  rho_ias = f_malloc((/ lr%d%n1i, lr%d%n2i, n3p, nmulti /),id='rho_ias')
  v_ias = f_malloc((/ lr%d%n1i, lr%d%n2i, n3p /),id='v_ias')

  if (nspin == 1) then
     ndipoles=2*nmulti
  else
     ndipoles=nmulti
  end if

  dipoles = f_malloc((/ 3, ndipoles /),id='dipoles')
  call to_zero(3*ndipoles,dipoles)

  !allocate coupling matrix elements
  K = f_malloc0((/ nmulti, nmulti /),id='K')

  !for nspin==1, define an auxiliary matrix for spin-off-diagonal terms (fxc part)
  if (nspin==1) then
     Kaux = f_malloc0((/ nmulti, nmulti /),id='Kaux')
  end if


  !exponent for Tamm-Dancoff Approximation
  if (tda) then
     ntda=0
  else
     ntda=1
  end if

  hfac=1.0_gp/(hxh*hyh*hzh)

  !we can start now to build the partial densities
  !and the corresponding partial potentials
  !bigger loop over a multiindex
  !the multindex is organised as I=a+(i-1)*norbv
  ik=0
  loop_i: do imulti=1,orbsvirt%norb*orbsocc%norb
     !calculate the orbital index
     iorbi=(imulti-1)/orbsvirt%norb+1
     iorba=imulti-(iorbi-1)*orbsvirt%norb

     !check the spin of the orbitals considered
     if (orbsocc%spinsgn(iorbi) == orbsvirt%spinsgn(iorba)) then
        ik=ik+1
        if (orbsocc%spinsgn(iorbi) == 1.0_gp) then
           ispin=1
        else if (orbsocc%spinsgn(iorbi) == -1.0_gp) then
           ispin=2
        end if
     else
        cycle loop_i
     end if

     do i3p=1,n3p
        ind3=(i3p-1)*lr%d%n1i*lr%d%n2i
        z=real(i3p+i3s-nbl3-1,gp)*hzh-chargec(3)
        do i2=1,lr%d%n2i
           ind2=(i2-1)*lr%d%n1i+ind3
           y=real(i2-nbl2-1,gp)*hyh-chargec(2)
           do i1=1,lr%d%n1i
              x=real(i1-nbl1-1,gp)*hxh-chargec(1)
              indi=i1+ind2+(iorbi-1)*lr%d%n1i*lr%d%n2i*n3p
              inda=i1+ind2+(iorba-1)*lr%d%n1i*lr%d%n2i*n3p
              rho_ias(i1,i2,i3p,ik)=hfac*psivirtr(inda)*psirocc(indi)
              !calculate dipole moments
              dipoles(1,ik)=dipoles(1,ik)+x*rho_ias(i1,i2,i3p,ik)
              dipoles(2,ik)=dipoles(2,ik)+y*rho_ias(i1,i2,i3p,ik)
              dipoles(3,ik)=dipoles(3,ik)+z*rho_ias(i1,i2,i3p,ik)
           end do
        end do
     end do
     if (.not. onlyfxc) then
        !copy the partial density onto the partial potential space to pass it to PSolver
        call vcopy(lr%d%n1i*lr%d%n2i*n3p,rho_ias(1,1,1,ik),1,v_ias(1,1,1),1)
        !partial potential term for each partial density
!        if (iproc == 0 .and. verbose > 1) then
!           write(*,*)'Poisson Solver application: orbitals (virt,occ):',iorba,iorbi
!        end if
        call H_potential('D',pkernel,v_ias(1,1,1),rho_ias,ehart,0.0_dp,.false.,&
             quiet='YES')
!        if (iproc ==0) write(*,*) 'ehart',ehart*2.0_gp
     !after the Poisson Solver we can calculate the Upper triangular part of the Coupling matrix
     end if

     jk=0
     loop_j: do jmulti=1,imulti
        !calculate the orbital index
        jorbi=(jmulti-1)/orbsvirt%norb+1
        jorba=jmulti-(jorbi-1)*orbsvirt%norb

        !check the spin of the orbitals considered
        if (orbsocc%spinsgn(jorbi) == orbsvirt%spinsgn(jorba)) then
           jk=jk+1
           if (orbsocc%spinsgn(jorbi) == 1.0_gp) then
              jspin=1
           else if (orbsocc%spinsgn(jorbi) == -1.0_gp) then
              jspin=2
           end if
        else
           cycle loop_j
        end if

        if (.not. onlyfxc) then
           !multiplication of the RPA part
           K(ik,jk)=hxh*hyh*hzh*&
                dot(lr%d%n1i*lr%d%n2i*n3p,rho_ias(1,1,1,jk),1,v_ias(1,1,1),1)
           !in the non spin-pol case the RPA part is the same for spin off diagonal
           if (nspin ==1) then
              Kaux(ik,jk)=K(ik,jk)
           end if

        end if
        !add the XC contribution

        !map the spin couples to the index of dvxcdrho, in the abinit convention
        if (dofxc) then
           index=ispin+jspin-1

           do i3p=1,n3p
              do i2=1,lr%d%n2i
                 do i1=1,lr%d%n1i
                    K(ik,jk)=K(ik,jk)+hxh*hyh*hzh*&
                         rho_ias(i1,i2,i3p,ik)*rho_ias(i1,i2,i3p,jk)*&
                         dvxcdrho(i1,i2,i3p,index)
                 end do
              end do
           end do

           !calculate the spin off-diagonal term for nspin==1
           if (nspin ==1) then
              index=2
              do i3p=1,n3p
                 do i2=1,lr%d%n2i
                    do i1=1,lr%d%n1i
                       Kaux(ik,jk)=Kaux(ik,jk)+hxh*hyh*hzh*&
                            rho_ias(i1,i2,i3p,ik)*rho_ias(i1,i2,i3p,jk)*&
                            dvxcdrho(i1,i2,i3p,index)
                    end do
                 end do
              end do
           end if


        end if
        !add factors from energy occupation numbers (for non-tda case)
        K(ik,jk)=K(ik,jk)*(2.0_wp*sqrt(orbsvirt%eval(iorba)-orbsocc%eval(iorbi))*&
             sqrt(orbsvirt%eval(jorba)-orbsocc%eval(jorbi)))**ntda
        if (nspin ==1) then
           Kaux(ik,jk)=Kaux(ik,jk)*(2.0_wp*sqrt(orbsvirt%eval(iorba)-orbsocc%eval(iorbi))*&
                sqrt(orbsvirt%eval(jorba)-orbsocc%eval(jorbi)))**ntda
        end if

     end do loop_j
  end do loop_i

  if (nproc > 1) then
     call mpiallred(K(1,1),nmulti**2,MPI_SUM,bigdft_mpi%mpi_comm)
     if (nspin ==1) call mpiallred(Kaux(1,1),nmulti**2,MPI_SUM,bigdft_mpi%mpi_comm)
     call mpiallred(dipoles(1,1),3*nmulti,MPI_SUM,bigdft_mpi%mpi_comm)
  end if

  if (nspin==1) then
     !copy the values of the dipoles in the second part of the array
     call vcopy(3*nmulti,dipoles(1,1),1,dipoles(1,nmulti+1),1)
     call dscal(3*ndipoles,hxh*hyh*hzh,dipoles(1,1),1)
  end if

  !equal the lower triangular part of the matrix
  do imulti=1,nmulti
     do jmulti=imulti+1,nmulti
        K(imulti,jmulti)=K(jmulti,imulti)
     end do
  end do

  if (nspin==1) then
     do imulti=1,nmulti
        do jmulti=imulti+1,nmulti
           Kaux(imulti,jmulti)=Kaux(jmulti,imulti)
        end do
     end do
  end if

  !add the A matrix to the diagonal part
  !overwrite K (Tamm-Dancoff Approx.)
  ik=0
  loop_i2: do imulti=1,orbsvirt%norb*orbsocc%norb
     !calculate the orbital index
     iorbi=(imulti-1)/orbsvirt%norb+1
     iorba=imulti-(iorbi-1)*orbsvirt%norb
     !check the spin of the orbitals considered
     if (orbsocc%spinsgn(iorbi) == orbsvirt%spinsgn(iorba)) then
        ik=ik+1
     else
        cycle loop_i2
     end if
     !energy difference of the eigenvalues
     K(ik,ik)=K(ik,ik)+(orbsvirt%eval(iorba)-orbsocc%eval(iorbi))**(ntda+1)
  end do loop_i2


  !print out the result
  ik=0
  loop_i3: do imulti=1,orbsvirt%norb*orbsocc%norb
     !calculate the orbital index
     iorbi=(imulti-1)/orbsvirt%norb+1
     iorba=imulti-(iorbi-1)*orbsvirt%norb
     if (orbsocc%spinsgn(iorbi) == orbsvirt%spinsgn(iorba)) then
        ik=ik+1
     else
        cycle loop_i3
     end if
!     if (iproc == 0) write(*,*) 'iorba,iorbi,K',iorba,iorbi,K(ik,ik)
  end do loop_i3

  call f_free(v_ias)


  if (tda) then

     !oscillator strength
     fi = f_malloc((/ 3, ndipoles /),id='fi')


     if (nspin == 1) then
        Kbig = f_malloc0((/ 2*nmulti, 2*nmulti /),id='Kbig')

        do ik=1,nmulti
           do jk=1,nmulti
              Kbig(ik,jk)=K(ik,jk)
              Kbig(ik+nmulti,jk+nmulti)=K(ik,jk)
              Kbig(ik+nmulti,jk)=Kaux(ik,jk)
              Kbig(ik,jk+nmulti)=Kaux(ik,jk)
           end do
        end do
        omega = f_malloc(2*nmulti,id='omega')

        lwork=6*nmulti !safe value
        work = f_malloc(lwork,id='work')


        call DSYEV('V','U',2*nmulti,Kbig,2*nmulti,omega,work,lwork,info)
        if (info /= 0) then
           call yaml_warning('Error, DSYEV' // trim(yaml_toa(info)))
           !write(*,*) 'Error, DSYEV',info
        end if

        !!treat the eigenvectors such that they have always positive numbers in the first elements
        !do jmulti=1,2*nmulti
        !   !if it is negative change the sign to all the values
        !   if (Kbig(1,jmulti) < 0.0_wp) then
        !      do imulti=1,2*nmulti
        !         Kbig(imulti,jmulti)=-Kbig(imulti,jmulti)
        !      end do
        !   end if
        !end do

        !transition dipoles
        call gemm('N','N',3,ndipoles,ndipoles,1.0_wp,dipoles(1,1),3,&
             Kbig(1,1),ndipoles,0.0_wp,fi(1,1),3)
        
        ! summary of the results and pretty printing
        if (iproc == 0) then
           if (tda) call yaml_comment('TAMM-DANCOFF APPROXIMATION',hfill='-')
           call yaml_open_sequence('Excitation Energy and Oscillator Strength')
           !write(6,'(/)')
           !if (tda) print *, 'TAMM-DANCOFF APPROXIMATION'
           !write(6,10)
!10         format ('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
           !write(6,20)
!20         format(t2,4x,'Excitation Energy',3x,'Oscillator Strength')
!           write(6,10)

           do imulti = 1, 2*nmulti
              call yaml_sequence(trim(yaml_toa((/ Ha_eV*omega(imulti),&
                    omega(imulti)*(2.0_gp/3.0_gp)*(fi(1,imulti)**2+fi(2,imulti)**2+fi(3,imulti)**2) /),&
                    fmt='(1pe10.3)')),advance='no')
              call yaml_comment(trim(yaml_toa(imulti,fmt='(i4.4)')))
              !write(6,30) imulti, Ha_eV*omega(imulti),omega(imulti)*(2./3.)*(fi(1,imulti)**2+fi(2,imulti)**2+fi(3,imulti)**2)
!30            format(t2,i3,2x,f9.4,12x,1pe10.3) 
           end do
           call yaml_close_sequence()

!          Extracting the excitation energies and Oscillator strength to plot absorption spectra
           open(unit=9, file='td_spectra.txt')
           write(9,'(a4)')'2  #(results in eV)' 
           do imulti = 1, min(100,2*nmulti) 
              write(9,'(f9.4,5x,1pe10.3)') Ha_eV*omega(imulti),&
                   omega(imulti)*(2.0_gp/3.0_gp)*(fi(1,imulti)**2+fi(2,imulti)**2+fi(3,imulti)**2)
           end do
           close(unit=9)

           !write(6,10)

           call yaml_open_sequence('Transition energies (eV)')
           do imulti = 1,2*nmulti
              call yaml_sequence(advance='no')
              call yaml_open_sequence(advance='no',flow=.true.)
              call yaml_map('Energy',trim(yaml_toa(Ha_eV*omega(imulti),fmt='(f10.5)')))
              !write(6,40)
!40            format('================================================')
              !write(6,50) imulti, Ha_eV*omega(imulti) 
!50            format(i3,2x,'Transition energy',2x,f10.5,1x,'eV')  
              !write(6,70)
!70            format('------------------------------------------------')
              ik=0
              do iorbi = 1, orbsocc%norb
                 do iorba = 1, orbsvirt%norb
                    jmulti =  (iorbi-1)*orbsvirt%norb+ iorba
                    if (abs(Kbig(jmulti,imulti)) > 5.d-02) then
                       if (ik /= 0) call yaml_newline()
                       ik = ik + 1
                       call yaml_open_map(flow=.true.)
                          call yaml_map('Transition',trim(yaml_toa((/ iorbi, iorba /))))
                          call yaml_map('Coeff',trim(yaml_toa(abs(Kbig(jmulti,imulti)),fmt='(1pe10.3)')))
                       call yaml_close_map()   
                       !write(6,60) iorbi, iorba,  abs(Kbig(jmulti,imulti))
!60                     format (i4,'----->',i3,2x,' Coeff =',1pe10.3) 
                    end if
                 end do
              end do

              !write(6,70)

              !ik = 0
              do iorbi = 1, orbsocc%norb
                 do iorba = 1, orbsvirt%norb
                    jmulti =  (iorbi-1)*orbsvirt%norb+ iorba
                    if (abs(Kbig(jmulti+nmulti,imulti)) > 5.d-02) then
                       if (ik /= 0) call yaml_newline()
                       ik = ik + 1
                       call yaml_open_map(flow=.true.)
                          call yaml_map('Transition',trim(yaml_toa((/ iorbi, iorba /))))
                          call yaml_map('Coeff',trim(yaml_toa(abs(Kbig(jmulti+nmulti,imulti)),fmt='(1pe10.3)')))
                       call yaml_close_map()
                       !write(6,60) iorbi, iorba,  abs(Kbig(jmulti+nmulti,imulti))
                    end if
                 end do
             end do
             call yaml_close_sequence(advance='no')
             call yaml_comment(trim(yaml_toa(imulti,fmt='(i4.4)')))

           end do
           call yaml_close_sequence()
           !write(6,40)

        end if

        call f_free(Kbig)
     
        call f_free(omega)
        omega = f_malloc(nmulti,id='omega')

        call f_free(work)
        lwork=3*nmulti !safe value
        work = f_malloc(lwork,id='work')


        !this second part is not needed
!!$        call DSYEV('V','U',nmulti,K,nmulti,omega,work,lwork,info)
!!$        if (info /= 0) then
!!$           call yaml_warning('Error, DSYEV' // trim(yaml_toa(info)))
!!$           !write(*,*) 'error, DSYEV',info
!!$        end if
!!$
!!$        !transition dipoles
!!$        call gemm('N','N',3,ndipoles,ndipoles,1.0_wp,dipoles(1,1),3,&
!!$             K(1,1),ndipoles,0.0_wp,fi(1,1),3)
!!$
!!$        !print eigenvalues
!!$        if (iproc == 0) then
!!$           call yaml_open_sequence('Excitation energies (Ha, eV, dipoles)')
!!$           do imulti=1,nmulti
!!$              call yaml_sequence( trim(yaml_toa( &
!!$                 & (/ omega(imulti),omega(imulti)*Ha_eV,fi(1,imulti),fi(2,imulti),fi(3,imulti) /),fmt='(f10.5)')),&
!!$                 & advance='no')
!!$              call yaml_comment(trim(yaml_toa(imulti,fmt='(i4.4)')))
!!$              !print '(a,i6,2(f10.5),3(f10.5))','excitation energies: Ha, eV, dipoles:' , &
!!$              !imulti,omega(imulti),omega(imulti)*Ha_eV,fi(1,imulti),fi(2,imulti),fi(3,imulti)
!!$           end do
!!$           call yaml_close_sequence()
!!$        end if
         
     end if
     
     call f_free(omega)
     call f_free(work)
  end if  

  call f_free(fi)
  call f_free(rho_ias)

  if (nspin ==1 ) then
     call f_free(Kaux)
  end if

  call f_free(K)
  call f_free(dipoles)

END SUBROUTINE coupling_matrix_prelim
