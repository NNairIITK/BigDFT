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
