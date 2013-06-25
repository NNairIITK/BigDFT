program wocc
   nup_tot=60
   ndw_tot=60

   nup_occ=35+4
   ndw_occ=35-4
   nup_vrt=nup_tot-nup_occ
   ndw_vrt=ndw_tot-ndw_occ

   write(*,*) nup_occ+nup_vrt, ndw_occ+ndw_vrt
   ii=0

   do i=1,nup_occ   
   ii=ii+1
   write(*,*) ii,1
   enddo

   do i=1,nup_vrt   
   ii=ii+1
   write(*,*) ii,0
   enddo

   do i=1,ndw_occ   
   ii=ii+1
   write(*,*) ii,1
   enddo

   do i=1,ndw_vrt   
   ii=ii+1
   write(*,*) ii,0
   enddo

end program wocc
