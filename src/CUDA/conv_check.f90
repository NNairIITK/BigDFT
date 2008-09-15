program conv_check
  use module_base
  implicit none
  integer  :: n1,n2,n3
  real(gp) :: hx,hy,hz,r2,sigma2,x,y,z,maxdiff,epot,arg
  real(wp), dimension(:,:,:), allocatable :: pot,psir,psi_in,psi_out
  !local variables
  character(len=*), parameter :: subname='conv_check'
  integer :: i,i_stat,i_all,j,i1,i2,i3,ntimes,ndat
  real(wp) :: tt
  real(gp) :: v,p,CPUtime,GPUtime
  real(kind=4), dimension(:,:,:), allocatable :: psi_cuda,v_cuda !temporary in view of wp 
  real(kind=4) :: t0,t1
  real(kind=8) :: psi_GPU,v_GPU,work_GPU !pointer to the GPU  memory addresses (with norb=1)
  integer, parameter :: lowfil1=-8,lupfil1=7 !for GPU computation
  integer, parameter :: lowfil2=-7,lupfil2=8 !for GPU computation
  real(kind=4) filCUDA1(lowfil1:lupfil1) !array of filters to be passed to CUDA interface
  data filCUDA1 / &
       8.4334247333529341094733325815816e-7_4,&
       -0.1290557201342060969516786758559028e-4_4,&
       0.8762984476210559564689161894116397e-4_4,&
       -0.30158038132690463167163703826169879e-3_4,&
       0.174723713672993903449447812749852942e-2_4,&
       -0.942047030201080385922711540948195075e-2_4,&
       0.2373821463724942397566389712597274535e-1_4,&
       0.612625895831207982195380597e-1_4,&
       0.9940415697834003993178616713_4,&
       -0.604895289196983516002834636e-1_4, &
       -0.2103025160930381434955489412839065067e-1_4,&
       0.1337263414854794752733423467013220997e-1_4,&
       -0.344128144493493857280881509686821861e-2_4,&
       0.49443227688689919192282259476750972e-3_4,&
       -0.5185986881173432922848639136911487e-4_4,&
       2.72734492911979659657715313017228e-6_4 /
  real(kind=4) filCUDA2(lowfil2:lupfil2) !array of filters to be passed to CUDA interface
  data filCUDA2 / &
       2.72734492911979659657715313017228e-6_4,&
       -0.5185986881173432922848639136911487e-4_4,&
       0.49443227688689919192282259476750972e-3_4,&
       -0.344128144493493857280881509686821861e-2_4,&
       0.1337263414854794752733423467013220997e-1_4,&
       -0.2103025160930381434955489412839065067e-1_4,&
       -0.604895289196983516002834636e-1_4,&
       0.9940415697834003993178616713_4,&
       0.612625895831207982195380597e-1_4,&
       0.2373821463724942397566389712597274535e-1_4,&
       -0.942047030201080385922711540948195075e-2_4,&
       0.174723713672993903449447812749852942e-2_4,&
       -0.30158038132690463167163703826169879e-3_4,&
       0.8762984476210559564689161894116397e-4_4,&
       -0.1290557201342060969516786758559028e-4_4,&
       8.4334247333529341094733325815816e-7_4 /

  !values for the grid
  n1=50
  n2=50
  n3=50
  ntimes=100
  hx=0.1e0_gp
  hy=0.1e0_gp
  hz=0.1e0_gp


  !set of one-dimensional convolutions
  !allocate arrays
  ndat=20000!(n2+1)*(n3+1)
  allocate(psi_in((2*n1+2),ndat,1+ndebug),stat=i_stat)
  call memocc(i_stat,psi_in,'psi_in',subname)
  allocate(psi_out(ndat,(2*n1+2),1+ndebug),stat=i_stat)
  call memocc(i_stat,psi_out,'psi_out',subname)

  !initialise array
  sigma2=0.25d0*(((2*n1+2)*hx)**2)
  do i=1,ndat
     do i1=1,2*n1+2
        x=hx*real(i1-n1-1,kind=8)
        !tt=abs(dsin(real(i1+i2+i3,kind=8)+.7d0))
        r2=x**2
        arg=0.5d0*r2/sigma2
        tt=dexp(-arg)

        psi_in(i1,i,1)=tt
     end do
  end do
 
  !take timings
  call cpu_time(t0)
  do i=1,ntimes
     call convrot_n_per(2*n1+1,ndat,psi_in,psi_out)
  end do
  call cpu_time(t1)

  CPUtime=real(t1-t0,kind=8)/real(ntimes,kind=8)


  allocate(psi_cuda((2*n1+2),ndat,1+ndebug),stat=i_stat)
  call memocc(i_stat,psi_cuda,'psi_cuda',subname)
  allocate(v_cuda(ndat,(2*n1+2),1+ndebug),stat=i_stat)
  call memocc(i_stat,v_cuda,'v_cuda',subname)


  !the input and output arrays must be reverted in this implementation
  v_cuda=real(psi_in,kind=4)
  do i=1,ndat
     do i1=1,2*n1+2
        v_cuda(i,i1,1)=real(psi_in(i1,i,1),kind=4)
        !write(16,'(2(i6),2(1pe24.17)')i,i1,v_cuda(i,i1,1),psi_in(i1,i,1)
     end do
  end do

  call GPU_allocate((2*n1+2)*ndat,psi_GPU,i_stat)
  call GPU_allocate((2*n1+2)*ndat,work_GPU,i_stat)

  !call GPU_send((2*n1+2)*ndat,v_cuda,psi_GPU,i_stat)
  call GPU_send((2*n1+2)*ndat,v_cuda,work_GPU,i_stat)

  !now the CUDA part
  !take timings
  call cpu_time(t0)
  do i=1,ntimes
     !call cuda C interface
     !call previous1dconv(2*n1+1,ndat,v_cuda,psi_cuda,filCUDA1,lowfil1,lupfil1)
     call m1dconv(2*n1+1,ndat,work_GPU,psi_GPU,filCUDA1,lowfil1,lupfil1)

  end do
  call cpu_time(t1)

  call GPU_receive((2*n1+2)*ndat,psi_cuda,psi_GPU,i_stat)

  call GPU_deallocate(psi_GPU,i_stat)
  call GPU_deallocate(work_GPU,i_stat)


  GPUtime=real(t1-t0,kind=8)/real(ntimes,kind=8)
  
  !check the differences between the results
  maxdiff=0.d0
  do i=1,ndat
     do i1=1,2*n1+2
        !write(17,'(2(i6),2(1pe24.17)')i,i1,psi_out(i,i1,1),psi_cuda(i1,i,1)
        maxdiff=max(abs(psi_out(i,i1,1)-real(psi_cuda(i1,i,1),kind=8)),maxdiff)
     end do
  end do

  print *,'timings,difference',CPUtime,GPUtime,maxdiff

  print *,'CPU/GPU ratio',CPUtime/GPUtime

  i_all=-product(shape(psi_in))
  deallocate(psi_in,stat=i_stat)
  call memocc(i_stat,i_all,'psi_in',subname)
  i_all=-product(shape(psi_out))
  deallocate(psi_out,stat=i_stat)
  call memocc(i_stat,i_all,'psi_out',subname)

  i_all=-product(shape(psi_cuda))
  deallocate(psi_cuda,stat=i_stat)
  call memocc(i_stat,i_all,'psi_cuda',subname)
  i_all=-product(shape(v_cuda))
  deallocate(v_cuda,stat=i_stat)
  call memocc(i_stat,i_all,'v_cuda',subname)



  !-------------------------3d case---------------------------
  ntimes=1
  !allocate arrays
  allocate(psi_in((2*n1+2),(2*n2+2),(2*n3+2)+ndebug),stat=i_stat)
  call memocc(i_stat,psi_in,'psi_in',subname)
  allocate(psi_out((2*n1+2),(2*n2+2),(2*n3+2)+ndebug),stat=i_stat)
  call memocc(i_stat,psi_out,'psi_out',subname)
  allocate(psir((2*n1+2),(2*n2+2),(2*n3+2)+ndebug),stat=i_stat)
  call memocc(i_stat,psir,'psir',subname)
  allocate(pot((2*n1+2),(2*n2+2),(2*n3+2)+ndebug),stat=i_stat)
  call memocc(i_stat,pot,'pot',subname)

  ! Wavefunction expressed everywhere in fine scaling functions 

  !fake initialisation, random numbers
  !here the grid spacings are the small ones
  sigma2=0.25d0*(((2*n1+2)*hx)**2+((2*n2+2)*hy)**2+((2*n3+2)*hz)**2)
  do i3=1,2*n3+2
     z=hz*real(i3-n3-1,kind=8)
     do i2=1,2*n2+2
        y=hy*real(i2-n2-1,kind=8)
        do i1=1,2*n1+2
           x=hx*real(i1-n1-1,kind=8)
           !tt=abs(dsin(real(i1+i2+i3,kind=8)+.7d0))
           r2=x**2+y**2+z**2
           arg=0.5d0*r2/sigma2
           tt=dexp(-arg)
           !same initialisation for psi and pot
           psi_in(i1,i2,i3)=tt
           pot(i1,i2,i3)=tt
        end do
     end do
  end do

  !take timings
  call cpu_time(t0)
  
  do j=1,ntimes
     !traditional convolution, CPU
     ! psir serves as a work array	   
     call convolut_magic_n_per(2*n1+1,2*n2+1,2*n3+1,psi_in,psir,psi_out) 

     ! Initialisation of potential energy  
     epot=0.0_gp
     do i3=1,2*n3+2
        do i2=1,2*n2+2
           do i1=1,2*n1+2
              v=real(pot(i1,i2,i3),gp)
              p=real(psir(i1,i2,i3),gp)
              tt=pot(i1,i2,i3)*psir(i1,i2,i3)
              epot=epot+p*v*p
              psir(i1,i2,i3)=tt
           end do
        end do
     end do

     call convolut_magic_t_per_self(2*n1+1,2*n2+1,2*n3+1,psir,psi_out)

  end do

  call cpu_time(t1)

  CPUtime=real(t1-t0,kind=8)/real(ntimes,kind=8)

  print *,'epot,CPUtime',epot,CPUtime

  !save the values of psi_out and of epot for GPU comparison

  !create the CUDA values
  allocate(psi_cuda((2*n1+2),(2*n2+2),(2*n3+2)+ndebug),stat=i_stat)
  call memocc(i_stat,psi_cuda,'psi_cuda',subname)
  allocate(v_cuda((2*n1+2),(2*n2+2),(2*n3+2)+ndebug),stat=i_stat)
  call memocc(i_stat,v_cuda,'v_cuda',subname)

  psi_cuda=real(psi_in,kind=4)
  v_cuda=real(pot,kind=4)
  
  ! Initialisation of potential energy  
  epot=0.0_gp

  !allocate the GPU memory
  call GPU_allocate((2*n1+2)*(2*n2+2)*(2*n3+2),psi_GPU,i_stat)
  call GPU_allocate((2*n1+2)*(2*n2+2)*(2*n3+2),v_GPU,i_stat)
  call GPU_allocate((2*n1+2)*(2*n2+2)*(2*n3+2),work_GPU,i_stat)

  call GPU_send((2*n1+2)*(2*n2+2)*(2*n3+2),psi_cuda,psi_GPU,i_stat)
  call GPU_send((2*n1+2)*(2*n2+2)*(2*n3+2),v_cuda,v_GPU,i_stat)

  !copy the data on GPU(must be separated)
  !call CUDA_ALLOC_MEM(1,2*n1+1,2*n2+1,2*n3+1,psi_cuda,v_cuda,psi_GPU,v_GPU,work_GPU)

  call cpu_time(t0)
  do i=1,ntimes

     !calculate the potential application on GPU
     call cuda_psi_to_vpsi(1,2*n1+1,2*n2+1,2*n3+1,psi_GPU,v_GPU,work_GPU,&
          filCUDA1,filCUDA2,lowfil1,lupfil1,lowfil2,lupfil2)

     !copy vpsi on the CPU
     call GPU_receive((2*n1+2)*(2*n2+2)*(2*n3+2),psi_cuda,psi_GPU,i_stat)
!     call cuda_fetch_vpsi(1,2*n1+1,2*n2+1,2*n3+1,psi_GPU,psi_cuda)

  end do
  call cpu_time(t1)

  GPUtime=real(t1-t0,kind=8)/real(ntimes,kind=8)

  !deallocate GPU memory
  call GPU_deallocate(psi_GPU,i_stat)
  call GPU_deallocate(v_GPU,i_stat)
  call GPU_deallocate(work_GPU,i_stat)

  !call CUDA_DEALLOCATE_MEM(1,psi_GPU,v_GPU,work_GPU)

  !check the differences between the results
  maxdiff=0.d0
  do i3=1,2*n3+2
     do i2=1,2*n2+2
        do i1=1,2*n1+2
           maxdiff=max(abs(psi_out(i1,i2,i3)-real(psi_cuda(i1,i2,i3),kind=8)),maxdiff)
        end do
     end do
  end do
 
  print *,'timings,difference',CPUtime,GPUtime,maxdiff
  
  i_all=-product(shape(psi_cuda))
  deallocate(psi_cuda,stat=i_stat)
  call memocc(i_stat,i_all,'psi_cuda',subname)
  i_all=-product(shape(v_cuda))
  deallocate(v_cuda,stat=i_stat)
  call memocc(i_stat,i_all,'v_cuda',subname)

  i_all=-product(shape(psi_in))
  deallocate(psi_in,stat=i_stat)
  call memocc(i_stat,i_all,'psi_in',subname)
  i_all=-product(shape(psi_out))
  deallocate(psi_out,stat=i_stat)
  call memocc(i_stat,i_all,'psi_out',subname)
  i_all=-product(shape(psir))
  deallocate(psir,stat=i_stat)
  call memocc(i_stat,i_all,'psir',subname)
  i_all=-product(shape(pot))
  deallocate(pot,stat=i_stat)
  call memocc(i_stat,i_all,'pot',subname)


end program conv_check
