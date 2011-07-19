
             program PSPimport
             !trivial conversion from old psp.par format to psppar pspcod 10
             implicit none
             real(8):: hsep(6,4,2),havg(6),hso(6),r_l(4),gpot(4),
     :                 rloc,rcov,rprb,rij,znuc,zion
             integer:: ng,l,lpx,lpj,i,j,nspol,ngpot
             character(22)::label1,label2
             character(1)::il(4)
             logical:: exists

             inquire(file='psp.par',exist=exists)
             if(.not.exists)then
                write(6,*)'This program is used to convert psp.par'
                write(6,*)'files to pseudopotentials in the ABINIT'
                write(6,*)'format (pspcod=10). No such file found.'
                stop
             end if

c                 read the data, code taken from previous version of pseudo

             write(6,*) '***        Reading data from psp.par      ***'
             open(unit=23,file='psp.par',form='formatted',
     1            status='unknown')
             read(23,*) ng,rij
             read(23,*) rcov,rprb 
             read(23,'(a)') label1
             read(23,'(a)') label2
             read(23,*) znuc,zion,rloc,gpot(1),gpot(2),gpot(3),gpot(4)
             read(23,*) lpx
             do l=0,lpx
                read(23,*) r_l(l+1),  (hsep(i,l+1,1),i=1,6)
                if (l.gt.0) read(23,*)(hsep(i,l+1,2),i=1,6)
             enddo
             close(23)



c            write out the psppar, code is nearly the same as in pseudo2.5
             nspol=1
             il(1)='s'
             il(2)='p'
             il(3)='d'
             il(4)='f'

             write(6,*)'The following output is appended to psppar ...'
             open(13,file='psppar',position='append')
             write(13,'(a)',advance='no')trim(label1)
             write( 6,'(a)',advance='no')trim(label1)
             write(13,'(i3,3(1x,g9.3),a)')ng,rij,rcov,rprb,
     :                            'ispp, ng rij rcov and rprb'
             write(06,'(i3,3(1x,g9.3),a)')ng,rij,rcov,rprb,
     :                            'ispp, ng rij rcov and rprb'

             write(13,'(2i5,2x,a)')int(znuc+.1),int(zion+.1),
     :                        ' zatom, zion'
             write( 6,'(2i5,2x,a)')int(znuc+.1),int(zion+.1),
     :                        ' zatom, zion'
             write(13,'(3x,i3,a,i2,a)')10,trim(label2),lpx,
     :                  ' 0 2002 0     pspcod,IXC,lmax,lloc,mmax,r2well'
             write( 6,'(3x,i3,a,i2,a)')10,trim(label2),lpx,
     :                ' 0 2002 0     pspcod,!IXC!,lmax,lloc,mmax,r2well'
             ngpot=0
             do j=1,4
                if (gpot(j).ne.0.d0) ngpot=j
             enddo
             if(ngpot==0)then
               write(13,'(e18.8,a)')rloc,' 0 rloc nloc ck (none)'
               write( 6,'(e18.8,a)')rloc,' 0 rloc nloc ck (none)'
             else
               write(label1,'(a,i1.1,a)')'(e16.8,i2,',ngpot,'e16.8,a)'
               write(13,label1)rloc,ngpot,gpot(1:ngpot),
     :                ' rloc nloc c1 .. cnloc'
               write( 6,label1)rloc,ngpot,gpot(1:ngpot),
     :                ' rloc nloc c1 .. cnloc'
             end if
             write(13,*)lpx+1, 'nnonloc'
             write( 6,*)lpx+1, 'nnonloc'
             do l=0,lpx
                do i=1,6
                   if (l.gt.1-nspol) then
                      havg(i)=((l+1)*hsep(i,l+1,1)+l*hsep(i,l+1,2))
     :                     /(2*l+1)
                      hso(i)=2*(hsep(i,l+1,1)-hsep(i,l+1,2))
     :                     /(2*l+1)
                   else
                      havg(i)=hsep(i,l+1,1)
                   endif
                enddo
c               get the matrix dimension lpj = 1,2 or 3
                lpj=1
                if(max(abs(havg(3)),abs(hso(3)))>1d-8)lpj=2
                if(max(abs(havg(6)),abs(hso(6)))>1d-8)lpj=3
                if(lpj==1)then
                    write(13,'(2x,e16.8,i3,e16.8,6x,a)')
     :                         r_l(l+1),lpj, havg(1)
     :                        ,il(l+1)//'-projector'
                    write( 6,'(2x,e16.8,i3,e16.8,6x,a)')
     :                         r_l(l+1),lpj, havg(1)
     :                        ,il(l+1)//'-projector'
                    if (l.gt.1-nspol)write(13,'(21x,e16.8)')hso(1)
                    if (l.gt.1-nspol)write( 6,'(21x,e16.8)')hso(1)
                elseif(lpj==2)then
                    write(13,'(2x,e16.8,i3,2e16.8,6x,a)')
     :                         r_l(l+1),lpj, havg(1:2)
     :                        ,il(l+1)//'-projector'
                    write( 6,'(2x,e16.8,i3,2e16.8,6x,a)')
     :                         r_l(l+1),lpj, havg(1:2)
     :                        ,il(l+1)//'-projector'
                    write(13,'(37x,e16.8)') havg(3)
                    write( 6,'(37x,e16.8)') havg(3)
                    if (l.gt.1-nspol)then
                      write(13,'(21x,2e16.8)') hso(1:2)
                      write( 6,'(21x,2e16.8)') hso(1:2)
                      write(13,'(37x, e16.8)') hso(3)
                      write( 6,'(37x, e16.8)') hso(3)
                    end if
                elseif(lpj==3)then
                    write(13,'(2x,e16.8,i3,3e16.8,6x,a)')
     :                         r_l(l+1),lpj,havg(1:2), havg(4)
     :                        ,il(l+1)//'-projector'
                    write( 6,'(2x,e16.8,i3,3e16.8,6x,a)')
     :                         r_l(l+1),lpj,havg(1:2), havg(4)
     :                        ,il(l+1)//'-projector'
                    write(13,'(37x,2e16.8)') havg(3),havg(5)
                    write( 6,'(37x,2e16.8)') havg(3),havg(5)
                    write(13,'(53x,e16.8)') havg(6)
                    write( 6,'(53x,e16.8)') havg(6)
                    if (l.gt.1-nspol)then
                      write(13,'(21x,3e16.8)') hso(1:2),hso(4)
                      write( 6,'(21x,3e16.8)') hso(1:2),hso(4)
                      write(13,'(37x,2e16.8)') hso(3),hso(5)
                      write( 6,'(37x,2e16.8)') hso(3),hso(5)
                      write(13,'(53x, e16.8)') hso(6)
                      write( 6,'(53x, e16.8)') hso(6)
                    end if
                end if
             enddo
             end program
