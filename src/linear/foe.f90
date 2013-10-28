!> @file
!! Fermi Operator Expansion Medthod
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Could still do more tidying - assuming all sparse matrices except for Fermi have the same pattern
subroutine foe(iproc, nproc, orbs, foe_obj, tmprtr, mode, &
           ham, ovrlp, fermi, ebs, itout, it_scc)
  use module_base
  use module_types
  use module_interfaces, except_this_one => foe
  use yaml_output
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nproc,itout,it_scc
  type(orbitals_data),intent(in) :: orbs
  type(foe_data),intent(inout) :: foe_obj
  real(kind=8),intent(inout) :: tmprtr
  integer,intent(in) :: mode
  type(sparseMatrix),intent(in) :: ovrlp, ham
  type(sparseMatrix),intent(inout) :: fermi
  real(kind=8),intent(out) :: ebs

  ! Local variables
  integer :: npl, istat, iall, jorb, info, ipl, i, it, ierr, ii, iiorb, jjorb, iseg, it_solver!, iorb
  integer :: isegstart, isegend, iismall, iseglarge, isegsmall, is, ie, iilarge, nsize_polynomial
  integer,parameter :: nplx=5000
  real(kind=8),dimension(:,:),allocatable :: cc, fermip, chebyshev_polynomials
  real(kind=8),dimension(:,:,:),allocatable :: penalty_ev
  real(kind=8) :: anoise, scale_factor, shift_value, sumn, sumnder, charge_diff, ef_interpol
  real(kind=8) :: evlow_old, evhigh_old, m, b, det, determinant, sumn_old, ef_old, bound_low, bound_up
  logical :: restart, adjust_lower_bound, adjust_upper_bound, calculate_SHS, interpolation_possible
  character(len=*),parameter :: subname='foe'
  real(kind=8),dimension(2) :: efarr, sumnarr, allredarr
  real(kind=8),dimension(:),allocatable :: hamscal_compr, ovrlpeff_compr, SHS
  real(kind=8),dimension(4,4) :: interpol_matrix, tmp_matrix
  real(kind=8),dimension(4) :: interpol_vector, interpol_solution
  integer,dimension(4) :: ipiv
  real(kind=8),parameter :: charge_tolerance=1.d-6 ! exit criterion

  !!real(8),dimension(100000) ::  work, eval, hamtmp, ovrlptmp


  call timing(iproc, 'FOE_auxiliary ', 'ON')

  ! initialization
  interpol_solution = 0.d0


  !!hamtmp(1:orbs%norb**2)=ham%matrix_compr
  !!ovrlptmp(1:orbs%norb**2)=ovrlp%matrix_compr
  !!call dsygv(1, 'v', 'l', orbs%norb, hamtmp, orbs%norb, ovrlptmp, orbs%norb, eval, work, 1000, ii)
  !!if (iproc==0) write(*,*) 'evals',eval(1), eval(orbs%norb)


  allocate(penalty_ev(orbs%norb,orbs%norbp,2), stat=istat)
  call memocc(istat, penalty_ev, 'penalty_ev', subname)


  allocate(ovrlpeff_compr(ovrlp%nvctr), stat=istat)
  call memocc(istat, ovrlpeff_compr, 'ovrlpeff_compr', subname)

  allocate(fermip(orbs%norb,orbs%norbp), stat=istat)
  call memocc(istat, fermip, 'fermip', subname)

  allocate(SHS(ovrlp%nvctr), stat=istat)
  call memocc(istat, SHS, 'SHS', subname)

  ii=0
  do iseg=1,ovrlp%nseg
      do jorb=ovrlp%keyg(1,iseg),ovrlp%keyg(2,iseg)
          iiorb = (jorb-1)/orbs%norb + 1
          jjorb = jorb - (iiorb-1)*orbs%norb
          ii=ii+1
          if (iiorb==jjorb) then
              ovrlpeff_compr(ii)=1.5d0-.5d0*ovrlp%matrix_compr(ii)
          else
              ovrlpeff_compr(ii)=-.5d0*ovrlp%matrix_compr(ii)
          end if
      end do  
  end do

  allocate(hamscal_compr(ham%nvctr), stat=istat)
  call memocc(istat, hamscal_compr, 'hamscal_compr', subname)




  evlow_old=1.d100
  evhigh_old=-1.d100

  if (mode==1) then
      
      stop 'not guaranteed to work anymore'

  else if (mode==2) then

      ! Size of one Chebyshev polynomial matrix in compressed form (distributed)
      nsize_polynomial=0
      if (orbs%norbp>0) then
          isegstart=fermi%istsegline(orbs%isorb_par(iproc)+1)
          if (orbs%isorb+orbs%norbp<orbs%norb) then
              isegend=fermi%istsegline(orbs%isorb_par(iproc+1)+1)-1
          else
              isegend=fermi%nseg
          end if
          !!$omp parallel default(private) shared(isegstart, isegend, fermi, nsize_polynomial)
          !!$omp do reduction(+:nsize_polynomial)
          do iseg=isegstart,isegend
              do jorb=fermi%keyg(1,iseg),fermi%keyg(2,iseg)
                  nsize_polynomial=nsize_polynomial+1
              end do
          end do
          !!$omp end do
          !!$omp end parallel
      end if


      ! Fake allocation, will be modified later
      allocate(chebyshev_polynomials(nsize_polynomial,1),stat=istat)
      call memocc(istat,chebyshev_polynomials,'chebyshev_polynomials',subname)

      ! Don't let this value become too small.
      foe_obj%bisection_shift = max(foe_obj%bisection_shift,1.d-4)

      efarr(1)=foe_obj%ef-foe_obj%bisection_shift
      efarr(2)=foe_obj%ef+foe_obj%bisection_shift
      sumnarr(1)=0.d0
      sumnarr(2)=1.d100


      adjust_lower_bound=.true.
      adjust_upper_bound=.true.

      calculate_SHS=.true.

      if (orbs%norbp>0) then
          call to_zero(orbs%norb*orbs%norbp, fermip(1,1))
      end if

      if (iproc==0) then
          call yaml_open_sequence('FOE to determine density kernel',label=&
               'it_foe'//trim(adjustl(yaml_toa(itout,fmt='(i3.3)')))//'-'//&
               trim(adjustl(yaml_toa(it_scc,fmt='(i3.3)'))))
      end if



      it=0
      it_solver=0
      main_loop: do 
          
          it=it+1

          if (iproc==0) then
              call yaml_newline()
              call yaml_sequence(advance='no')
              call yaml_open_map(flow=.true.)
              call yaml_comment('it FOE:'//yaml_toa(it,fmt='(i6)'),hfill='-')
          end if
          
          if (adjust_lower_bound) then
              foe_obj%ef=efarr(1)
          else if (adjust_upper_bound) then
              foe_obj%ef=efarr(2)
          end if
      

          ! Scale the Hamiltonian such that all eigenvalues are in the intervall [-1:1]
          if (foe_obj%evlow/=evlow_old .or. foe_obj%evhigh/=evhigh_old) then
              scale_factor=2.d0/(foe_obj%evhigh-foe_obj%evlow)
              shift_value=.5d0*(foe_obj%evhigh+foe_obj%evlow)
              ii=0
              do iseg=1,ovrlp%nseg
                  do jorb=ovrlp%keyg(1,iseg),ovrlp%keyg(2,iseg)
                      ii=ii+1
                      hamscal_compr(ii)=scale_factor*(ham%matrix_compr(ii)-shift_value*ovrlp%matrix_compr(ii))
                  end do  
              end do
              calculate_SHS=.true.
          else
              calculate_SHS=.false.
          end if
          evlow_old=foe_obj%evlow
          evhigh_old=foe_obj%evhigh

          !!foe_obj%ef = foe_obj%evlow+1.d-4*it

          ! Determine the degree of the polynomial
          npl=nint(3.0d0*(foe_obj%evhigh-foe_obj%evlow)/foe_obj%fscale)
          if (npl>nplx) stop 'npl>nplx'

          ! Array the holds the Chebyshev polynomials. Needs to be recalculated
          ! every thime the Hamiltonian has been modified.
          if (calculate_SHS) then
              iall=-product(shape(chebyshev_polynomials))*kind(chebyshev_polynomials)
              deallocate(chebyshev_polynomials,stat=istat)
              call memocc(istat,iall,'chebyshev_polynomials',subname)
              allocate(chebyshev_polynomials(nsize_polynomial,npl),stat=istat)
              call memocc(istat,chebyshev_polynomials,'chebyshev_polynomials',subname)
          end if

          if (iproc==0) then
              !!write( *,'(1x,a,i0)') repeat('-',75 - int(log(real(it))/log(10.))) // ' FOE it=', it
              !!write(*,'(1x,a,2x,i0,2es12.3,es18.9,3x,i0)') 'FOE: it, evlow, evhigh, efermi, npl', &
              !!                                              it, foe_obj%evlow, foe_obj%evhigh, foe_obj%ef, npl
              !!write(*,'(1x,a,2x,2es13.5)') 'Bisection bounds: ', efarr(1), efarr(2)
              call yaml_map('lower bisection bound',efarr(1),fmt='(es10.3)')
              call yaml_map('upper bisection bound',efarr(1),fmt='(es10.3)')
              call yaml_map('lower eigenvalue bound',foe_obj%evlow,fmt='(es10.3)')
              call yaml_map('upper eigenvalue bound',foe_obj%evhigh,fmt='(es10.3)')
              call yaml_map('pol degree',npl)
              call yaml_map('eF',foe_obj%ef,fmt='(es19.12)')
              call yaml_map('conv crit',charge_tolerance,fmt='(es10.3)')
          end if


          allocate(cc(npl,3), stat=istat)
          call memocc(istat, cc, 'cc', subname)

          if (foe_obj%evlow>=0.d0) then
              stop 'ERROR: lowest eigenvalue must be negative'
          end if
          if (foe_obj%evhigh<=0.d0) then
              stop 'ERROR: highest eigenvalue must be positive'
          end if

          call timing(iproc, 'FOE_auxiliary ', 'OF')
          call timing(iproc, 'chebyshev_coef', 'ON')

          call chebft(foe_obj%evlow, foe_obj%evhigh, npl, cc(1,1), foe_obj%ef, foe_obj%fscale, tmprtr)
          call chder(foe_obj%evlow, foe_obj%evhigh, cc(1,1), cc(1,2), npl)
          call chebft2(foe_obj%evlow, foe_obj%evhigh, npl, cc(1,3))
          call evnoise(npl, cc(1,3), foe_obj%evlow, foe_obj%evhigh, anoise)

          call timing(iproc, 'chebyshev_coef', 'OF')
          call timing(iproc, 'FOE_auxiliary ', 'ON')
        
          !!if (iproc==0) then
          !!    call pltwght(npl,cc(1,1),cc(1,2),foe_obj%evlow,foe_obj%evhigh,foe_obj%ef,foe_obj%fscale,tmprtr)
          !!    call pltexp(anoise,npl,cc(1,3),foe_obj%evlow,foe_obj%evhigh)
          !!end if
        
        
          if (orbs%nspin==1) then
              do ipl=1,npl
                  cc(ipl,1)=2.d0*cc(ipl,1)
                  cc(ipl,2)=2.d0*cc(ipl,2)
                  cc(ipl,3)=2.d0*cc(ipl,3)
              end do
          end if
        
        
          call timing(iproc, 'FOE_auxiliary ', 'OF')

          if (calculate_SHS) then
              ! sending it ovrlp just for sparsity pattern, still more cleaning could be done
              !!if (iproc==0) write(*,*) 'Need to recalculate the Chebyshev polynomials.'
              if (iproc==0) call yaml_map('Chebyshev polynomials','recalculated')
              call chebyshev_clean(iproc, nproc, npl, cc, orbs, foe_obj, ovrlp, fermi, hamscal_compr, &
                   ovrlpeff_compr, calculate_SHS, nsize_polynomial, SHS, fermip, penalty_ev, chebyshev_polynomials)
          else
              ! The Chebyshev polynomials are already available
              !!if (iproc==0) write(*,*) 'Can use the Chebyshev polynomials from memory.'
              if (iproc==0) call yaml_map('Chebyshev polynomials','from memory')
              call chebyshev_fast(iproc, nsize_polynomial, npl, orbs, fermi, chebyshev_polynomials, cc, fermip)
          end if 

          call timing(iproc, 'FOE_auxiliary ', 'ON')


          restart=.false.

          ! Check the eigenvalue bounds. Only necessary if calculate_SHS is true
          ! (otherwise this has already been checked in the previous iteration).
          if (calculate_SHS) then
              ! The penalty function must be smaller than the noise.
              bound_low=0.d0
              bound_up=0.d0
              if (orbs%norbp>0) then
                  isegstart=ovrlp%istsegline(orbs%isorb_par(iproc)+1)
                  if (orbs%isorb+orbs%norbp<orbs%norb) then
                      isegend=ovrlp%istsegline(orbs%isorb_par(iproc+1)+1)-1
                  else
                      isegend=ovrlp%nseg
                  end if
                  !$omp parallel default(private) shared(isegstart, isegend, orbs, penalty_ev, ovrlp, bound_low, bound_up)
                  !$omp do reduction(+:bound_low,bound_up)
                  do iseg=isegstart,isegend
                      ii=ovrlp%keyv(iseg)-1
                      do jorb=ovrlp%keyg(1,iseg),ovrlp%keyg(2,iseg)
                          ii=ii+1
                          iiorb = (jorb-1)/orbs%norb + 1
                          jjorb = jorb - (iiorb-1)*orbs%norb
                          bound_low = bound_low + penalty_ev(jjorb,iiorb-orbs%isorb,2)*ovrlp%matrix_compr(ii)
                          bound_up = bound_up +penalty_ev(jjorb,iiorb-orbs%isorb,1)*ovrlp%matrix_compr(ii)
                      end do  
                  end do
                  !$omp end do
                  !$omp end parallel
              end if

              allredarr(1)=bound_low
              allredarr(2)=bound_up
              call mpiallred(allredarr, 2, mpi_sum, bigdft_mpi%mpi_comm, ierr)
              if (allredarr(1)>anoise) then
                  if (iproc==0) then
                      !!write(*,'(1x,a,2es12.3)') 'WARNING: lowest eigenvalue to high; penalty function, noise: ', &
                      !!                          allredarr(1), anoise
                      !!write(*,'(1x,a)') 'Increase magnitude by 20% and cycle.'
                      call yaml_map('lowest eigenvalue within bounds','No; increase magnitude by 20% and cycle')
                  end if
                  foe_obj%evlow=foe_obj%evlow*1.2d0
                  restart=.true.
              else
                  if (iproc==0) then
                      !!write(*,'(1x,a)') 'Lowest eigenvalue within interval, can continue.'
                      call yaml_map('lowest eigenvalue within bounds','Yes; can continue')
                  end if
              end if
              if (allredarr(2)>anoise) then
                  if (iproc==0) then
                      !!write(*,'(1x,a,2es12.3)') 'WARNING: highest eigenvalue to low; penalty function, noise: ', &
                      !!                          allredarr(2), anoise
                      !!write(*,'(1x,a)') 'Increase magnitude by 20% and cycle.'
                      call yaml_map('highest eigenvalue within bounds','No; increase magnitude by 20% and cycle')
                  end if
                  foe_obj%evhigh=foe_obj%evhigh*1.2d0
                  restart=.true.
              else
                  if (iproc==0) then
                      !!write(*,'(1x,a)') 'Highest eigenvalue within interval, can continue.'
                      call yaml_map('highest eigenvalue within bounds','Yes; can continue')
                  end if
              end if
          end if

          iall=-product(shape(cc))*kind(cc)
          deallocate(cc, stat=istat)
          call memocc(istat, iall, 'cc', subname)

          if (restart) then
              if (iproc==0) then
                  call yaml_close_map()
                  call bigdft_utils_flush(unit=6)
              end if
              cycle
          end if
              
        

          sumn=0.d0
          sumnder=0.d0
          if (orbs%norbp>0) then
              isegstart=ovrlp%istsegline(orbs%isorb_par(iproc)+1)
              if (orbs%isorb+orbs%norbp<orbs%norb) then
                  isegend=ovrlp%istsegline(orbs%isorb_par(iproc+1)+1)-1
              else
                  isegend=ovrlp%nseg
              end if
              !!$omp parallel default(private) shared(isegstart, isegend, orbs, fermip, ovrlp, sumn) 
              !!$omp do reduction(+:sumn)
              do iseg=isegstart,isegend
                  ii=ovrlp%keyv(iseg)-1
                  do jorb=ovrlp%keyg(1,iseg),ovrlp%keyg(2,iseg)
                      ii=ii+1
                      iiorb = (jorb-1)/orbs%norb + 1
                      jjorb = jorb - (iiorb-1)*orbs%norb
                      sumn = sumn + fermip(jjorb,iiorb-orbs%isorb)*ovrlp%matrix_compr(ii)
                  end do  
              end do
              !!$omp end do
              !!$omp end parallel
          end if

          if (nproc>1) then
              call mpiallred(sumn, 1, mpi_sum, bigdft_mpi%mpi_comm, ierr)
          end if



          ! Make sure that the bounds for the bisection are negative and positive
          charge_diff = sumn-foe_obj%charge
          if (adjust_lower_bound) then
              if (iproc==0) call yaml_map('checking lower bisection bound, charge diff',charge_diff,fmt='(es9.2)')
              if (charge_diff<=0.d0) then
                  ! Lower bound okay
                  adjust_lower_bound=.false.
                  foe_obj%bisection_shift=foe_obj%bisection_shift*9.d-1
                  sumnarr(1)=sumn
                  !!it_solver=it_solver+1
                  !!ii=min(it_solver,4)
                  !!do i=1,4
                  !!    interpol_matrix(ii,1)=foe_obj%ef**3
                  !!    interpol_matrix(ii,2)=foe_obj%ef**2
                  !!    interpol_matrix(ii,3)=foe_obj%ef
                  !!    interpol_matrix(ii,4)=1
                  !!end do
                  !!interpol_vector(ii)=(sumn-foe_obj%charge)
                  !!if (iproc==0) write(*,'(1x,a)') 'lower bisection bound is okay.'
                  if (iproc==0) then
                      call yaml_map('lower bound okay',.true.)
                      call yaml_close_map()
                      call bigdft_utils_flush(unit=6)
                  end if
                  cycle !now check the upper bound
              else
                  efarr(1)=efarr(1)-foe_obj%bisection_shift
                  foe_obj%bisection_shift=foe_obj%bisection_shift*1.1d0
                  !!if (iproc==0) write(*,'(1x,a,es12.5)') &
                  !!    'lower bisection bound does not give negative charge difference: diff=',charge_diff
                  if (iproc==0) then
                      call yaml_map('lower bound okay',.false.)
                      call yaml_close_map()
                      call bigdft_utils_flush(unit=6)
                  end if
                  cycle !move the bisection bound
              end if
          end if
          if (adjust_upper_bound) then
              if (iproc==0) call yaml_map('checking upper bisection bound, charge diff',charge_diff,fmt='(es9.2)')
              if (charge_diff>=0.d0) then
                  ! Upper bound okay
                  adjust_upper_bound=.false.
                  foe_obj%bisection_shift=foe_obj%bisection_shift*9.d-1
                  sumnarr(2)=sumn
                  !!if (iproc==0) write(*,'(1x,a)') 'upper bisection bound is okay.'
                  if (iproc==0) call yaml_map('upper bound okay',.true.)
              else
                  efarr(2)=efarr(2)+foe_obj%bisection_shift
                  foe_obj%bisection_shift=foe_obj%bisection_shift*1.1d0
                  !!if (iproc==0) write(*,'(1x,a,es12.5)') &
                  !!    'upper bisection bound does not give positive charge difference: diff=',charge_diff
                  if (iproc==0) then
                      call yaml_map('upper bound okay',.false.)
                      call yaml_close_map()
                      call bigdft_utils_flush(unit=6)
                      cycle !move the bisection bound
                  end if
              end if
          end if

          !!! Cycle if one of the two bounds is not okay
          !!if (adjust_lower_bound .or. adjust_upper_bound) cycle

          it_solver=it_solver+1

          ! Check whether the system behaves reasonably.
          interpolation_possible=.true.
          if (it_solver>1) then
              if (iproc==0) call yaml_newline()
              if (iproc==0) call yaml_open_map('Interpolation check',flow=.true.)
              if (iproc==0) call yaml_map('D eF',foe_obj%ef-ef_old,fmt='(es13.6)')
              if (iproc==0) call yaml_map('D Tr',sumn-sumn_old,fmt='(es13.6)')
              if (foe_obj%ef>ef_old .and. sumn<sumn_old) then
                  !!if (iproc==0) then
                  !!    !!write(*,'(1x,a,2es13.5)') 'WARNING: Fermi energy was raised, but the trace still decreased: Deltas=',&
                  !!    !!                  foe_obj%ef-ef_old,sumn-sumn_old
                  !!    !!write(*,'(1x,a)') 'Cubic interpolation not possible under this circumstances.'
                  !!end if
                  interpolation_possible=.false.
              end if
              if (foe_obj%ef<ef_old .and. sumn>sumn_old) then
                  !!if (iproc==0) then
                  !!    !!write(*,'(1x,a,2es13.5)') 'WARNING: Fermi energy was lowered, but the trace still increased: Deltas=',&
                  !!    !!                  foe_obj%ef-ef_old,sumn-sumn_old
                  !!    !!write(*,'(1x,a)') 'Cubic interpolation not possible under this circumstances.'
                  !!end if
                  interpolation_possible=.false.
              end if
              if (foe_obj%ef>ef_old .and. sumn<sumn_old .or. foe_obj%ef<ef_old .and. sumn>sumn_old) then
                  if (iproc==0) call yaml_map('interpol possible',.false.)
              else
                  if (iproc==0) call yaml_map('interpol possible',.true.)
              end if
              if (iproc==0) call yaml_close_map()
              if (iproc==0) call yaml_newline()
          end if
          if (.not.interpolation_possible) then
              ! Set the history for the interpolation to zero.
              it_solver=0
          end if

          ef_old=foe_obj%ef
          sumn_old=sumn

          ! Shift up the old results.
          if (it_solver>4) then
              do i=1,4
                  interpol_matrix(1,i)=interpol_matrix(2,i)
                  interpol_matrix(2,i)=interpol_matrix(3,i)
                  interpol_matrix(3,i)=interpol_matrix(4,i)
              end do
              interpol_vector(1)=interpol_vector(2)
              interpol_vector(2)=interpol_vector(3)
              interpol_vector(3)=interpol_vector(4)
          end if
          !LG: if it_solver==0 this index comes out of bounds!
          ii=max(min(it_solver,4),1)
          interpol_matrix(ii,1)=foe_obj%ef**3
          interpol_matrix(ii,2)=foe_obj%ef**2
          interpol_matrix(ii,3)=foe_obj%ef
          interpol_matrix(ii,4)=1
          interpol_vector(ii)=sumn-foe_obj%charge

          ! Solve the linear system interpol_matrix*interpol_solution=interpol_vector
          if (it_solver>=4) then
              do i=1,ii
                  interpol_solution(i)=interpol_vector(i)
                  tmp_matrix(i,1)=interpol_matrix(i,1)
                  tmp_matrix(i,2)=interpol_matrix(i,2)
                  tmp_matrix(i,3)=interpol_matrix(i,3)
                  tmp_matrix(i,4)=interpol_matrix(i,4)
              end do
              !!if (iproc==0) then
              !!    do i=1,ii
              !!        write(*,'(4es14.6,5x,es14.6)') tmp_matrix(i,1:4), interpol_solution(i)
              !!    end do
              !!end if

              call dgesv(ii, 1, tmp_matrix, 4, ipiv, interpol_solution, 4, info)
              if (info/=0) then
                 if (iproc==0) write(*,'(1x,a,i0)') 'ERROR in dgesv (FOE), info=',info
              end if
              !!if (iproc==0) write(*,'(a,4es14.7)') 'interpol_solution(1:4)',interpol_solution(1:4)


              call get_roots_of_cubic_polynomial(interpol_solution(1), interpol_solution(2), &
                   interpol_solution(3), interpol_solution(4), foe_obj%ef, ef_interpol)
          end if


          ! Adjust the bounds for the bisection.
          if (charge_diff<0.d0) then
              efarr(1)=foe_obj%ef
              sumnarr(1)=sumn
          else if (charge_diff>=0.d0) then
              efarr(2)=foe_obj%ef
              sumnarr(2)=sumn
          end if


          ! Calculate the new Fermi energy.
          !!if (iproc==0) write(*,'(a,i6,2es16.6)') 'it_solver, abs(sumn-charge), foe_obj%ef_interpol_chargediff', it_solver, abs(sumn-foe_obj%charge), foe_obj%ef_interpol_chargediff
          !!if (iproc==0) write(*,'(a,5es16.6)') 'efarr(1), efarr(2), sumnarr(1), sumnarr(2), charge', efarr(1), efarr(2), sumnarr(1), sumnarr(2), foe_obj%charge
          if (iproc==0) call yaml_newline()
          if (iproc==0) call yaml_open_map('Search new eF',flow=.true.)
          if (it_solver>=4 .and.  abs(sumn-foe_obj%charge)<foe_obj%ef_interpol_chargediff) then
              det=determinant(iproc,4,interpol_matrix)
              !!if (iproc==0) write(*,'(1x,a,2es10.2)') 'determinant of interpolation matrix, limit:', &
              !!                                       det, foe_obj%ef_interpol_det
              if (iproc==0) call yaml_map('det',det,fmt='(es10.3)')
              if (iproc==0) call yaml_map('limit',foe_obj%ef_interpol_det,fmt='(es10.3)')
              if(abs(det)>foe_obj%ef_interpol_det) then
                  foe_obj%ef=ef_interpol
                  !!if (iproc==0) write(*,'(1x,a)') 'new fermi energy from cubic interpolation'
                  if (iproc==0) call yaml_map('method','cubic interpolation')
              else
                  ! linear interpolation
                  !!if (iproc==0) write(*,'(1x,a)') 'new fermi energy from linear interpolation'
                  if (iproc==0) call yaml_map('method','linear interpolation')
                  m = (interpol_vector(4)-interpol_vector(3))/(interpol_matrix(4,3)-interpol_matrix(3,3))
                  b = interpol_vector(4)-m*interpol_matrix(4,3)
                  foe_obj%ef = -b/m
              end if
          else
              ! Use mean value of bisection and secant method
              ! Secant method solution
              foe_obj%ef = efarr(2)-(sumnarr(2)-foe_obj%charge)*(efarr(2)-efarr(1))/(sumnarr(2)-sumnarr(1))
              ! Add bisection solution
              foe_obj%ef = foe_obj%ef + .5d0*(efarr(1)+efarr(2))
              ! Take the mean value
              foe_obj%ef=.5d0*foe_obj%ef
              !!if (iproc==0) write(*,'(1x,a)') 'new fermi energy from bisection / secant method'
              if (iproc==0) call yaml_map('method','bisection / secant method')
          end if
          if (iproc==0) then
              call yaml_close_map()
              call yaml_newline()
          end if


          if (iproc==0) then
              !!write(*,'(1x,a,2es21.13)') 'trace of the Fermi matrix, derivative matrix:', sumn, sumnder
              !!write(*,'(1x,a,2es13.4)') 'charge difference, exit criterion:', sumn-foe_obj%charge, charge_tolerance
              !!write(*,'(1x,a,es21.13)') 'suggested Fermi energy for next iteration:', foe_obj%ef
              call yaml_map('iter',it)
              call yaml_map('Tr(K)',sumn,fmt='(es16.9)')
              call yaml_map('charge diff',sumn-foe_obj%charge,fmt='(es16.9)')
          end if

          if (iproc==0) then
              call yaml_close_map()
              call bigdft_utils_flush(unit=6)
          end if

          if (abs(charge_diff)<charge_tolerance) then
              exit
          end if


      end do main_loop

      !!if (iproc==0) then
      !!    write( *,'(1x,a,i0)') repeat('-',84 - int(log(real(it))/log(10.)))
      !!end if

      if (iproc==0) then
          call yaml_close_sequence()
      end if




  end if


  call timing(iproc, 'FOE_auxiliary ', 'OF')
  call timing(iproc, 'chebyshev_comm', 'ON')

  call to_zero(fermi%nvctr, fermi%matrix_compr(1))

  if (orbs%norbp>0) then
      isegstart=fermi%istsegline(orbs%isorb_par(iproc)+1)
      if (orbs%isorb+orbs%norbp<orbs%norb) then
          isegend=fermi%istsegline(orbs%isorb_par(iproc+1)+1)-1
      else
          isegend=fermi%nseg
      end if
      !$omp parallel default(private) shared(isegstart, isegend, orbs, fermip, fermi)
      !$omp do
      do iseg=isegstart,isegend
          ii=fermi%keyv(iseg)-1
          do jorb=fermi%keyg(1,iseg),fermi%keyg(2,iseg)
              ii=ii+1
              iiorb = (jorb-1)/orbs%norb + 1
              jjorb = jorb - (iiorb-1)*orbs%norb
              fermi%matrix_compr(ii)=fermip(jjorb,iiorb-orbs%isorb)
          end do
      end do
      !$omp end do
      !$omp end parallel
  end if

  call mpiallred(fermi%matrix_compr(1), fermi%nvctr, mpi_sum, bigdft_mpi%mpi_comm, ierr)


  call timing(iproc, 'chebyshev_comm', 'OF')
  call timing(iproc, 'FOE_auxiliary ', 'ON')



  scale_factor=1.d0/scale_factor
  shift_value=-shift_value

  ebs=0.d0
  iismall=0
  iseglarge=1
  do isegsmall=1,ham%nseg
      do
          is=max(ham%keyg(1,isegsmall),fermi%keyg(1,iseglarge))
          ie=min(ham%keyg(2,isegsmall),fermi%keyg(2,iseglarge))
          iilarge=fermi%keyv(iseglarge)-fermi%keyg(1,iseglarge)
          do i=is,ie
              iismall=iismall+1
              ebs = ebs + fermi%matrix_compr(iilarge+i)*hamscal_compr(iismall)
          end do
          if (ie>=is) exit
          iseglarge=iseglarge+1
      end do
  end do
  ebs=ebs*scale_factor-shift_value*sumn


  iall=-product(shape(chebyshev_polynomials))*kind(chebyshev_polynomials)
  deallocate(chebyshev_polynomials,stat=istat)
  call memocc(istat,iall,'chebyshev_polynomials',subname)

  iall=-product(shape(penalty_ev))*kind(penalty_ev)
  deallocate(penalty_ev, stat=istat)
  call memocc(istat, iall, 'penalty_ev', subname)

  iall=-product(shape(ovrlpeff_compr))*kind(ovrlpeff_compr)
  deallocate(ovrlpeff_compr, stat=istat)
  call memocc(istat, iall, 'ovrlpeff_compr', subname)

  iall=-product(shape(hamscal_compr))*kind(hamscal_compr)
  deallocate(hamscal_compr, stat=istat)
  call memocc(istat, iall, 'hamscal_compr', subname)

  iall=-product(shape(fermip))*kind(fermip)
  deallocate(fermip, stat=istat)
  call memocc(istat, iall, 'fermip', subname)

  iall=-product(shape(SHS))*kind(SHS)
  deallocate(SHS, stat=istat)
  call memocc(istat, iall, 'SHS', subname)


  call timing(iproc, 'FOE_auxiliary ', 'OF')
     
      deallocate(penalty_ev, stat=istat)


end subroutine foe




! Calculates chebychev expansion of fermi distribution.
! Taken from numerical receipes: press et al
subroutine chebft(A,B,N,cc,ef,fscale,tmprtr)
  use module_base, pi => pi_param
  implicit none
  
  ! Calling arguments
  real(kind=8),intent(in) :: A, B, ef, fscale, tmprtr
  integer,intent(in) :: n
  real(8),dimension(n),intent(out) :: cc

  ! Local variables
  integer :: k, j
  real(kind=8) :: bma, bpa, y, arg, fac, tt, erfcc
  real(kind=8),dimension(5000) :: cf
  !real(kind=8),parameter :: pi=4.d0*atan(1.d0)

  if (n>5000) stop 'chebft'
  bma=0.5d0*(b-a)
  bpa=0.5d0*(b+a)
  do k=1,n
      y=cos(pi*(k-0.5d0)*(1.d0/n))
      arg=y*bma+bpa
      if (tmprtr.eq.0.d0) then
          cf(k)=.5d0*erfcc((arg-ef)*(1.d0/fscale))
      else
          cf(k)=1.d0/(1.d0+exp( (arg-ef)*(1.d0/tmprtr) ) )
      end if
  end do
  fac=2.d0/n
  do j=1,n
      tt=0.d0
      do  k=1,n
          tt=tt+cf(k)*cos((pi*(j-1))*((k-0.5d0)*(1.d0/n)))
      end do
      cc(j)=fac*tt
  end do

end subroutine chebft



! Calculates chebychev expansion of fermi distribution.
! Taken from numerical receipes: press et al
subroutine chebft2(a,b,n,cc)
  use module_base, pi => pi_param
  implicit none

  ! Calling arguments
  real(kind=8),intent(in) :: a, b
  integer,intent(in) :: n
  real(kind=8),dimension(n),intent(out) :: cc

  ! Local variables
  integer :: k, j
  !real(kind=8),parameter :: pi=4.d0*atan(1.d0)
  real(kind=8) :: tt, y, arg, fac, bma, bpa
  real(kind=8),dimension(5000) :: cf

  if (n>5000) stop 'chebft2'
  bma=0.5d0*(b-a)
  bpa=0.5d0*(b+a)
  ! 3 gives broder safety zone than 4
  !tt=3.0d0*n/(b-a)
  tt=4.d0*n/(b-a)
  do k=1,n
      y=cos(pi*(k-0.5d0)*(1.d0/n))
      arg=y*bma+bpa
      cf(k)=exp((arg-b)*tt)
  end do
  fac=2.d0/n
  do j=1,n
      tt=0.d0
      do k=1,n
          tt=tt+cf(k)*cos((pi*(j-1))*((k-0.5d0)*(1.d0/n)))
      end do
      cc(j)=fac*tt
  end do
end subroutine chebft2


! Calculates chebychev expansion of the derivative of Fermi distribution.
subroutine chder(a,b,c,cder,n)
  implicit none

  ! Calling arguments
  real(kind=8),intent(in) :: a, b
  integer,intent(in) :: n
  real(8),dimension(n),intent(in) :: c
  real(8),dimension(n),intent(out) :: cder

  ! Local variables
  integer :: j
  real(kind=8) :: con

  cder(n)=0.d0
  cder(n-1)=2*(n-1)*c(n)
  if(n>=3)then
      do j=n-2,1,-1
        cder(j)=cder(j+2)+2*j*c(j+1)
      end do
  end if
  con=2.d0/(b-a)
  do j=1,n
      cder(j)=cder(j)*con
  end do

end subroutine chder


!> Determine noise level
subroutine evnoise(npl,cc,evlow,evhigh,anoise)
  implicit none
  
  ! Calling arguments
  integer,intent(in) :: npl
  real(kind=8),dimension(npl),intent(in) :: cc
  real(kind=8),intent(in) :: evlow, evhigh
  real(kind=8),intent(out) :: anoise
  
  ! Local variables
  real(kind=8) :: fact, dist, ddx, cent, tt, x, chebev
  
  fact=1.d0
  dist=(fact*evhigh-fact*evlow)
  ddx=dist/(10*npl)
  cent=.5d0*(fact*evhigh+fact*evlow)
  !!tt=abs(chebev(evlow,evhigh,npl,cent,cc))
  !!do x=ddx,.25d0*dist,ddx
  !!    tt=max(tt,abs(chebev(evlow,evhigh,npl,cent+x,cc)), &
  !!       & abs(chebev(evlow,evhigh,npl,cent-x,cc)))
  !!end do
  ! Rewritten version ob the above loop
  tt=abs(chebev(evlow,evhigh,npl,cent,cc))
  x=ddx
  do 
      tt=max(tt,abs(chebev(evlow,evhigh,npl,cent+x,cc)), &
         & abs(chebev(evlow,evhigh,npl,cent-x,cc)))
      x=x+ddx
      if (x>=.25d0*dist) exit
  end do
  !anoise=1.d0*tt
  anoise=20.d0*tt

end subroutine evnoise



!> Calculates the error function complement with an error of less than 1.2E-7
function erfcc(x)
  implicit none

  ! Calling arguments
  real(8),intent(in) :: x
  real(8) :: erfcc

  ! Local variables
  real(8) :: z, t

  z=abs(x)
  t=1.d0/(1.+0.5d0*z)
  erfcc=t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+ &
        & t*(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+ &
        & t*(1.48851587+t*(-.82215223+t*.17087277)))))))))
  if (x.lt.0.) erfcc=2.D0-erfcc

end function erfcc


!> Evaluates chebychev expansion
function chebev(a,b,m,x,cc)
  implicit none
  
  ! Calling arguments
  real(kind=8),intent(in) :: a, b, x
  integer,intent(in) :: m
  real(kind=8),dimension(m),intent(in) :: cc
  real(kind=8) :: chebev
  
  ! Local variables
  integer :: j
  real(kind=8) :: d, dd, y, sv
  
  d=0.d0
  dd=0.d0
  y=2.d0*(2.d0*x-a-b)/(b-a)
  do j=m,2,-1
      sv=d
      d=y*d-dd+cc(j)
      dd=sv
  end do
  chebev= -dd + 0.5d0*(y*d+cc(1))

end function chebev
 



! plots the approximate fermi distribution
        subroutine pltwght(npl,cc,cder,evlow,evhigh,ef,fscale,tmprtr)
          implicit none

          ! Calling arguments
          integer,intent(in) :: npl
          real(kind=8),dimension(npl),intent(in) :: cc, cder
          real(kind=8),intent(in) :: evlow, evhigh, ef, fscale, tmprtr

          ! Local variables
          integer :: ic
          real(kind=8) :: ddx, x, tt, err, chebev

        open (unit=66,file='fermi',status='unknown')
!     header for my favourite plotting program
        write(66,*) ' 3'
        write(66,*) ' #LINETYPE{132}'
65        format(a,f5.2,a,i3,a)
        write(66,65) ' #TITLE{WEIGHT DISTR. for fscale=', fscale,' npl=',npl,'}'
        write(66,*) ' #XCAPT{ENERGY in eV}'
        write(66,*) ' #XYAPT{WEIGHT DISTR.}'
        write(66,*) ' #2YAXIS{2}'
        write(66,*) ' #YLOGSCALE2'
        write(66,*) ' #YCAPT2{ERROR}'
        write(66,*) ' $'
!
! plot chebechev expansion of weight distribution function
!
        ddx=(evhigh-evlow)/(10*npl)
! number of plot p[oints
        ic=0
        !!do x=evlow,evhigh,ddx
        !!    ic=ic+1
        !!end do
        x=evlow
        do
            ic=ic+1
            x=x+ddx
            if (x>=evhigh) exit
        end do
! weight distribution
        write(66,*) ic
        !!do x=evlow,evhigh,ddx
        !!    write(66,*) x,CHEBEV(evlow,evhigh,npl,x,cc)
        !!end do
        x=evlow
        do
            write(66,*) x,CHEBEV(evlow,evhigh,npl,x,cc)
            x=x+ddx
            if (x>=evhigh) exit
        end do
! derivative
        write(66,*) ic
        !!do x=evlow,evhigh,ddx
        !!    write(66,*) x,-CHEBEV(evlow,evhigh,npl,x,cder)
        !!end do
        x=evlow
        do
            write(66,*) x,-CHEBEV(evlow,evhigh,npl,x,cder)
            x=x+ddx
            if (x>=evhigh) exit
        end do
! error
        write(66,*) ic
        !!do x=evlow,evhigh,ddx
        !!    tt=tmprtr
        !!    if (tmprtr.eq.0.d0) tt=1.d-16
        !!    err=CHEBEV(evlow,evhigh,npl,x,cc) -1.d0/(1.d0+exp((x-ef)/tt))
        !!    write(66,*) x,err
        !!end do
        x=evlow
        do
            tt=tmprtr
            if (tmprtr.eq.0.d0) tt=1.d-16
            err=CHEBEV(evlow,evhigh,npl,x,cc) -1.d0/(1.d0+exp((x-ef)/tt))
            write(66,*) x,err
            x=x+ddx
            if (x>=evhigh) exit
        end do

        close(unit=66)
end subroutine pltwght




! plots the approximate fermi distribution
subroutine pltexp(anoise,npl,cc,evlow,evhigh)
        implicit none

        ! Calling arguments
        integer,intent(in) :: npl
        real(kind=8),dimension(npl),intent(in) :: cc
        real(kind=8),intent(in) :: anoise, evlow, evhigh

        ! Local variables
        integer :: ic
        real(kind=8) :: fact, ddx, tt, chebev, x

        open (unit=66,file='exp',status='unknown')
!     header for my favourite plotting program
        write(66,*) ' 2'
        write(66,*) ' #LINETYPE{12}'
        write(66,*) ' #TITLE{exponential}'
        write(66,*) ' #YLOGSCALE'
        write(66,*) ' #XCAPT{ENERGY in eV}'
        write(66,*) ' $'
!
        fact=1.25d0
! plot chebechev expansion of weight distribution function
!
        ddx=(fact*evhigh-fact*evlow)/(10*npl)
! number of plot p[oints
        ic=0
        !!do x=fact*evlow,fact*evhigh,ddx
        !!    ic=ic+1
        !!end do
        x=fact*evlow
        do
            ic=ic+1
            x=x+ddx
            if (x>=fact*evhigh) exit
        end do
! first curve
        write(66,*) ic
        !!do x=fact*evlow,fact*evhigh,ddx
        !!    tt=CHEBEV(evlow,evhigh,npl,x,cc)
        !!    if (abs(tt).lt.anoise) tt=anoise
        !!    write(66,*) x,tt
        !!end do
        x=fact*evlow
        do
            tt=CHEBEV(evlow,evhigh,npl,x,cc)
            if (abs(tt).lt.anoise) tt=anoise
            write(66,*) x,tt
            x=x+ddx
            if (x>=fact*evhigh) exit
        end do
! second curve
        write(66,*) ic
        !!do x=fact*evhigh,fact*evlow,-ddx
        !!    tt=CHEBEV(evlow,evhigh,npl,x,cc)
        !!    if (abs(tt).lt.anoise) tt=anoise
        !!    write(66,*) fact*evhigh-(x-fact*evlow),tt
        !!end do
        x=fact*evhigh
        do
            tt=CHEBEV(evlow,evhigh,npl,x,cc)
            if (abs(tt).lt.anoise) tt=anoise
            write(66,*) fact*evhigh-(x-fact*evlow),tt
            x=x-ddx
            if (x<=fact*evlow) exit
        end do

        close(unit=66)
end subroutine pltexp



! Finds the real root of the equation ax**3 + bx**2 + cx + d which is closest to target_solution
subroutine get_roots_of_cubic_polynomial(a, b, c, d, target_solution, solution)
  use module_base
  implicit none

  ! Calling arguments
  real(kind=8),intent(in) :: a, b, c, d
  real(kind=8),intent(in) :: target_solution
  real(kind=8),intent(out) :: solution

  ! Local variables
  complex(kind=8) :: a_c, b_c, c_c, d_c, Q_c, S_c, ttp_c, ttm_c
  complex(kind=8),dimension(3) :: sol_c
  real(kind=8) :: ttmin, tt
  integer :: i

  a_c=cmplx(a,0.d0,kind=8)
  b_c=cmplx(b,0.d0,kind=8)
  c_c=cmplx(c,0.d0,kind=8)
  d_c=cmplx(d,0.d0,kind=8)

  Q_c = sqrt( (2*b_c**3-9*a_c*b_c*c_c+27*a_c**2*d_c)**2 - 4*(b_c**2-3*a_c*c_c)**3 )
  S_c = ( .5d0*(Q_c+2*b_c**3-9*a_c*b_c*c_c+27*a_c**2*d_c) )**(1.d0/3.d0)
  ttp_c = cmplx(1.d0,sqrt(3.d0),kind=8)
  ttm_c = cmplx(1.d0,-sqrt(3.d0),kind=8)
  sol_c(1) = -b_c/(3*a_c) - S_c/(3*a_c) - (b_c**2-3*a_c*c_c)/(3*a_c*S_c)
  sol_c(2) = -b_c/(3*a_c) + (S_c*ttp_c)/(6*a_c) + ttm_c*(b_c**2-3*a_c*c_c)/(6*a_c*S_c)
  sol_c(3) = -b_c/(3*a_c) + (S_c*ttm_c)/(6*a_c) + ttp_c*(b_c**2-3*a_c*c_c)/(6*a_c*S_c)
  !!if (iproc==0) then
  !!    write(*,*) 'sol 1', sol_c(1)
  !!    write(*,*) 'sol 2', sol_c(2)
  !!    write(*,*) 'sol 3', sol_c(3)
  !!end if

  ! Select the real solution that is closest to target_solution
  ttmin=1.d100
  do i=1,3
      if (abs(aimag(sol_c(i)))>1.d-14) cycle !complex solution
      tt=abs(real(sol_c(i))-target_solution)
      if (tt<ttmin) then
          ttmin=tt
          solution=real(sol_c(i))
      end if
  end do

end subroutine get_roots_of_cubic_polynomial



real(kind=8) function determinant(iproc, n, mat)
    implicit none

    ! Calling arguments
    integer,intent(in) :: iproc, n
    real(kind=8),dimension(n,n),intent(in) :: mat

    ! Local variables
    integer :: i, info
    integer,dimension(n) :: ipiv
    real(kind=8),dimension(n,n) :: mat_tmp
    real(kind=8) :: sgn

    call dcopy(n**2, mat, 1, mat_tmp, 1)

    call dgetrf(n, n, mat_tmp, n, ipiv, info)
    if (info/=0) then
        if (iproc==0) write(*,'(a,i0,a)') 'ERROR in dgetrf, info=',info,'. Set determinant to zero.'
        determinant=0
        return
    end if

    determinant=1.d0
    do i=1,n
        determinant=determinant*mat_tmp(i,i)
    end do

    sgn=1.d0
    do i=1,n
        if(ipiv(i)/=i) then
            sgn=-sgn
        end if
    end do

    determinant=sgn*determinant   

end function determinant


subroutine compress_polynomial_vector(iproc, nsize_polynomial, orbs, fermi, vector, vector_compressed)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nsize_polynomial
  type(orbitals_data),intent(in) :: orbs
  type(sparseMatrix),intent(in) :: fermi
  real(kind=8),dimension(orbs%norb,orbs%norbp),intent(in) :: vector
  real(kind=8),dimension(nsize_polynomial),intent(out) :: vector_compressed

  ! Local variables
  integer :: isegstart, isegend, iseg, ii, jorb, iiorb, jjorb


  if (orbs%norbp>0) then
      isegstart=fermi%istsegline(orbs%isorb_par(iproc)+1)
      if (orbs%isorb+orbs%norbp<orbs%norb) then
          isegend=fermi%istsegline(orbs%isorb_par(iproc+1)+1)-1
      else
          isegend=fermi%nseg
      end if
      ii=0
      !!$omp parallel default(private) shared(isegstart, isegend, orbs, fermi, vector, vector_compressed)
      !!$omp do
      do iseg=isegstart,isegend
          !ii=fermi%keyv(iseg)-1
          do jorb=fermi%keyg(1,iseg),fermi%keyg(2,iseg)
              ii=ii+1
              iiorb = (jorb-1)/orbs%norb + 1
              jjorb = jorb - (iiorb-1)*orbs%norb
              vector_compressed(ii)=vector(jjorb,iiorb-orbs%isorb)
              !write(300,*) 'ii, jjorb, iiorb-orbs%isorb', ii, jjorb, iiorb-orbs%isorb
          end do
      end do
      !!$omp end do
      !!$omp end parallel
  end if
end subroutine compress_polynomial_vector



subroutine uncompress_polynomial_vector(iproc, nsize_polynomial, orbs, fermi, vector_compressed, vector)
  use module_base
  use module_types
  implicit none

  ! Calling arguments
  integer,intent(in) :: iproc, nsize_polynomial
  type(orbitals_data),intent(in) :: orbs
  type(sparseMatrix),intent(in) :: fermi
  real(kind=8),dimension(nsize_polynomial),intent(in) :: vector_compressed
  real(kind=8),dimension(orbs%norb,orbs%norbp),intent(out) :: vector

  ! Local variables
  integer :: isegstart, isegend, iseg, ii, jorb, iiorb, jjorb


  if (orbs%norbp>0) then
      call to_zero(orbs%norb*orbs%norbp, vector(1,1))
      isegstart=fermi%istsegline(orbs%isorb_par(iproc)+1)
      if (orbs%isorb+orbs%norbp<orbs%norb) then
          isegend=fermi%istsegline(orbs%isorb_par(iproc+1)+1)-1
      else
          isegend=fermi%nseg
      end if
      ii=0
      !!$omp parallel default(private) shared(isegstart, isegend, orbs, fermi, vector, vector_compressed)
      !!$omp do
      do iseg=isegstart,isegend
          !ii=fermi%keyv(iseg)-1
          do jorb=fermi%keyg(1,iseg),fermi%keyg(2,iseg)
              ii=ii+1
              iiorb = (jorb-1)/orbs%norb + 1
              jjorb = jorb - (iiorb-1)*orbs%norb
              vector(jjorb,iiorb-orbs%isorb)=vector_compressed(ii)
              !write(*,*) 'ii, iiorb-orbs%isorb, jjorb', ii, iiorb-orbs%isorb, jjorb
          end do
      end do
      !!$omp end do
      !!$omp end parallel
  end if
end subroutine uncompress_polynomial_vector
