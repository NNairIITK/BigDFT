
      subroutine pj2test(hsep,lpx,lpmx,lmx,nspin,nsmx,r_l,is)

      implicit none
      integer lpx,lpmx,ispin,nspin,nsmx,i,j,l,lmx,ll
      real*8 hsep(6,lpmx,nsmx),r_l(lmx),ohsep(6,lpmx,nsmx)
      real*8 h11,h12,h13,h22,h23,h33,hh11,hh12,hh13,hh22,hh23,hh33
      character*7 is(2)

c     save hsep()-values
      do ll=0,lpx
         do i=1,min(2*ll+1,nspin)
            do j=1,6
               ohsep(j,ll+1,i)= hsep(j,ll+1,i)
            enddo
         enddo
      enddo

c     give hij and projectors in nonorthogonal and othonormal space!

      write(6,*) 'test projectors....'
      write(6,*)'--------------------------------------------------'
      if (lpx.ge.0) then
         write(6,*) 'representation with nonorthogonal  projectors ',
     :        'p_i(l,r):'
         write(6,'(i4,t60,a)') lpx ,
     :        'lpx, (Projectors for l=0..lpx)'
         do l=0,lpx
            write(6,*) 'l=',l
            write(6,'(f7.3,t8,6e11.3,t76,a)') r_l(l+1),
     :           (hsep(i,l+1,1),i=1,6),'r_l(),hsep(), '//is(1)
            if (l.gt.0 .and. nspin.eq.2) 
     :           write(6,'(t8,6e11.3,t76,a)') 
     :           (hsep(i,l+1,2),i=1,6),'       hsep(), '//is(2)
            write(6,*)
            if (l.eq.0) then
               write(6,*)'p1(l=0,r)=2/(E**(r**2/(2*rl**2))*Pi**(1/4)',
     :              '*rl**(3/2))'
               write(6,*)'p2(l=0,r)=4*r**2/(Sqrt(15)*E**(r**2/(2*rl',
     :              '**2))*Pi**(1/4)*rl**(7/2))'
               write(6,*)'p3(l=0,r)=8*r**4/(3*Sqrt(105)*E**(r**2/',
     :              '(2*rl**2))*Pi**(1/4)*rl**(11/2))'
            else if (l.eq.1) then
               write(6,*) 'p1(l=1,r)=2*Sqrt(2/3)*r/(E**(r**2/',
     :              '(2*rl**2))*Pi**(1/4)*rl**(5/2))'
               write(6,*) 'p2(l=1,r)=4*Sqrt(2/105)*r**3/(E**(r**2/(2',
     :              '*rl**2))*Pi**(1/4)*rl**(9/2))'
               write(6,*) 'p3(l=1,r)=8*Sqrt(2/1155)*r**5/(3*E**(r**',
     :              '2/(2*rl**2))*Pi**(1/4)*rl**(13/2))'
            else if (l.eq.2) then
               write(6,*) 'p1(l=2,r)=4*r**2/(Sqrt(15)*E**(r**2/(2*rl',
     :              '**2))*Pi**(1/4)*rl**(7/2))'
               write(6,*) 'p2(l=2,r)=8*r**4/(3*Sqrt(105)*E**(r**2/',
     :              '(2*rl**2))*Pi**(1/4)*rl**(11/2))'
               write(6,*) 'p3(l=2,r)=16*r**6/(3*Sqrt(15015)*E**(r**2',
     :              '/(2*rl**2))*Pi**(1/4)*rl**(15/2))'
            else if (l.eq.3) then
               write(6,*) 'p1(l=3,r)=4*Sqrt(2/105)*r**3/(E**(r**2/',
     :              '(2*rl**2))*Pi**(1/4)*rl**(9/2))'
               write(6,*) 'p2(l=3,r)=8*Sqrt(2/1155)*r**5/(3*E**(r**2/',
     :              '(2*rl**2))*Pi**(1/4)*rl**(13/2))'
               write(6,*)'p3(l=3,r)=16*Sqrt(2/1001)*r**7/(45*E**',
     :              '(r**2/(2*rl**2))*Pi**(1/4)*rl**(17/2))'
            endif
            write(6,*)
         enddo


         write(6,*)'--------------------------------------------------'
         write(6,*) 'transformation to orthonormal  projectors P_i(l,r)'

         do ll=0,lpx
c            write(6,*)'ll=',ll
            do ispin=1,min(2*ll+1,nspin)
c               write(6,*)'ispin=',ispin
               h11=hsep(1,ll+1,ispin)
               h12=hsep(2,ll+1,ispin)
               h22=hsep(3,ll+1,ispin)
               h13=hsep(4,ll+1,ispin)
               h23=hsep(5,ll+1,ispin)
               h33=hsep(6,ll+1,ispin)
c               write(6,*)'hij=',h11,h12,h22
               if (ll.eq.0) then

      HH11=h11 + 1.549193338482967d0*h12 + 0.975900072948533d0*h13 
     :     + 0.6d0*h22 + 0.7559289460184545d0*h23 
     :     + 0.2380952380952381d0*h33
      HH12=0.6324555320336759d0*h12 + 0.7968190728895957d0*h13 + 
     :     0.4898979485566356d0*h22 + 0.925820099772551d0*h23 
     :     + 0.3888078956798695d0*h33
      HH13=0.3563483225498992d0*h13 + 0.2760262237369417d0*h23 + 
     :     0.173880176985767d0*h33
      HH22=0.4d0*h22 + 1.007905261357939d0*h23+0.6349206349206349d0*h33
      HH23=0.02839451399999733d0*(7.937253933193773d0*h23 + 10.d0*h33)
      HH33=0.126984126984127d0*h33

                 elseif (ll.eq.1) then

      HH11=h11 + 1.690308509457033d0*h12 + 1.189176780021126d0*h13
     :     + 0.7142857142857144d0*h22 + 1.005037815259212d0*h23 
     :     + 0.3535353535353535d0*h33
      HH12=0.5345224838248488d0*h12 + 0.7521014330903549d0*h13 
     :     + 0.4517539514526256d0*h22 + 0.953462589245592d0*h23 
     :     + 0.4471907802258314d0*h33
      HH13=0.2842676218074805d0*h13 + 0.240249990052149d0*h23 
     :     + 0.1690222275826415d0*h33
      HH22=0.2857142857142857d0*h22 + 0.80403025220737d0*h23 + 
     :     0.5656565656565657d0*h33
      HH23=0.01527129183875666d0*(9.9498743710662d0*h23 + 14.d0*h33)
      HH33=0.0808080808080808d0*h33

                 elseif (ll.eq.2) then

      HH11=h11 + 1.763834207376394d0*h12 + 1.327493036606129d0*h13 + 
     :     0.7777777777777778d0*h22 + 1.170738814009927d0*h23 
     :     + 0.4405594405594406d0*h33
      HH12=0.4714045207910317d0*h12 + 0.7095748751868991d0*h13 + 
     :     0.4157397096415491d0*h22 + 0.938679328162116d0*h23 + 
     :     0.4709778528806361d0*h33
      HH13=0.236524958395633d0*h13 + 0.2085954062582479d0*h23 + 
     :     0.1569926176268787d0*h33
      HH12=0.4714045207910317d0*h12 + 0.7095748751868991d0*h13 + 
     :     0.4157397096415491d0*h22 + 0.938679328162116d0*h23 + 
     :     0.4709778528806361d0*h33
      HH22=0.2222222222222222d0*h22 + 0.6689936080056727d0*h23 + 
     :     0.5034965034965035d0*h33
      HH23=0.00932400932400932d0*(11.9582607431014d0*h23 + 18.d0*h33)
      HH33=0.05594405594405595d0*h33

                 elseif (ll.eq.3) then

      HH11=h11 + 1.809068067466582d0*h12 + 1.425050606388851d0*h13 + 
     :     0.818181818181818d0*h22 + 1.289006773270979d0*h23 + 
     :     0.5076923076923077d0*h33
      HH12=0.0006593070220853591d0*(646.741834119303d0*h12 + 
     :     1018.911183568028d0*h13 + 585.d0*h22 + 
     :     1382.459764333125d0*h23 + 726.d0*h33)
      HH13=0.2025478734167333d0*h13 + 0.1832114449657378d0*h23 + 
     :     0.144320484917644d0*h33
      HH22=0.1818181818181818d0*h22 + 0.5728918992315464d0*h23 + 
     :     0.4512820512820513d0*h33
      HH23=0.006184848093902844d0*(13.96424004376894d0*h23+22.d0*h33)
      HH33=0.04102564102564103d0*h33

c               write(6,*)'HHij=',hh11,hh12,hh22
               endif
c               write(6,*)'hsep1',hsep(1,ll+1,ispin),
c     :              hsep(2,ll+1,ispin),hsep(3,ll+1,ispin)
               hsep(1,ll+1,ispin)=HH11
               hsep(2,ll+1,ispin)=HH12
               hsep(3,ll+1,ispin)=HH22
               hsep(4,ll+1,ispin)=HH13
               hsep(5,ll+1,ispin)=HH23
               hsep(6,ll+1,ispin)=HH33
c               write(6,*)'hsep1',hsep(1,ll+1,ispin),
c     :              hsep(2,ll+1,ispin),hsep(3,ll+1,ispin)
            enddo
         enddo

         write(6,'(i4,t60,a)') lpx ,
     :        'lpx, (Projectors for l=0..lpx)'
         do l=0,lpx
            write(6,*) 'l=',l
            write(6,'(f7.3,t8,6e11.3,t76,a)') r_l(l+1),
     :           (hsep(i,l+1,1),i=1,6),'r_l(),hsep(), '//is(1)
            if (l.gt.0 .and. nspin.eq.2) 
     :           write(6,'(t8,6e11.3,t76,a)') 
     :           (hsep(i,l+1,2),i=1,6),'       hsep(), '//is(2)
            write(6,*)
            if (l.eq.0) then
               write(6,*)'P1(l=0,r)=',
     :              '                  p1(l=0,r)'
               write(6,*)'P2(l=0,r)=',
     :              ' -Sqrt(3/2)     * p1(l=0,r)',
     :              ' +Sqrt(5/2)     * p2(l=0,r)'
               write(6,*)'P3(l=0,r)=',
     :              ' Sqrt(15/2)/2   * p1(l=0,r)',
     :              ' -5/Sqrt(2)     * p2(l=0,r)',
     :              ' +3*Sqrt(7/2)/2 * p3(l=0,r)'
            else if (l.eq.1) then
               write(6,*)'P1(l=1,r)=',
     :              '                    p1(l=1,r)'
               write(6,*)'P2(l=1,r)=',
     :              ' -Sqrt(5/2)       * p1(l=1,r)',
     :              ' +Sqrt(7/2)       * p2(l=1,r)'
               write(6,*)'P3(l=1,r)=',
     :              ' Sqrt(35/2)/2    * p1(l=1,r)',
     :              ' -7/Sqrt(2)      * p2(l=1,r)',
     :              ' +3*Sqrt(11/2)/2 * p3(l=1,r)'

            else if (l.eq.2) then
               write(6,*)'P1(l=2,r)=',
     :              '                  p1(l=2,r)'
               write(6,*)'P2(l=2,r)=',
     :              ' -Sqrt(7/2)     * p1(l=2,r)',
     :              ' +3/Sqrt(2)     * p2(l=2,r)'
               write(6,*)'P3(l=2,r)=',
     :              ' 3*Sqrt(7/2)/2  * p1(l=2,r)',
     :              ' -9/Sqrt(2)     * p2(l=2,r)',
     :              ' +Sqrt(143/2)/2 * p3(l=2,r)'

            else if (l.eq.3) then
               write(6,*)'P1(l=3,r)=',
     :              '                  p1(l=3,r)'
               write(6,*)'P2(l=3,r)=',
     :              ' -3/Sqrt(2)     * p1(l=3,r)',
     :              ' +Sqrt(11/2)    * p2(l=3,r)'
               write(6,*)'P3(l=3,r)=',
     :              ' 3*Sqrt(11/2)/2 * p1(l=3,r)',
     :              ' -11/Sqrt(2)    * p2(l=3,r)',
     :              ' +Sqrt(195/2)/2 * p3(l=3,r)'
            endif
            write(6,*)
         enddo

         write(6,*)'---------------------------------------------------'
         write(6,*) 'back transformation assuming orthogonal projectors'
         write(6,*) 'differences to original values:'

         do ll=0,lpx
            do ispin=1,min(2*ll+1,nspin)
               h11=hsep(1,ll+1,ispin)
               h22=hsep(3,ll+1,ispin)
               h33=hsep(6,ll+1,ispin)
               
               if (ll.eq.0) then

               HH11=H11 + 1.5d0*H22 + 1.875d0*H33
               HH12=-1.936491673103709d0*H22 - 4.841229182759272d0*H33
               HH13=3.842606537234849d0*H33
               HH22=2.5d0*H22 + 12.5d0*H33
               HH23=-9.92156741649221d0*H33
               HH33=7.875d0*H33
            

               elseif (ll.eq.1) then

               HH11=H11 + 2.5d0*H22 + 4.375d0*H33
               HH12=-2.958039891549808d0*H22 - 10.35313962042433d0*H33
               HH13=7.358031326380719d0*H33
               HH22=3.5d0*H22 + 24.5d0*H33
               HH23=-17.41228014936585d0*H33
               HH33=12.375d0*H33

               elseif (ll.eq.2) then

               HH11=H11 + 3.5d0*H22 + 7.875d0*H33
               HH12=-3.968626966596886d0*H22-17.85882134968598d0*H33
               HH13=11.86446901466728d0*H33
               HH22=4.5d0*H22 + 40.5d0*H33
               HH23=-26.90608667197814d0*H33
               HH33=17.875d0*H33

               elseif (ll.eq.3) then

               HH11=H11 + 4.5d0*H22 + 12.375d0*H33
               HH12=-4.9749371855331d0*H22 - 27.36215452043205d0*H33
               HH13=17.36780426536412d0*H33
               HH22=5.5d0*H22 + 60.5d0*H33
               HH23=-38.40166012036458d0*H33
               HH33=24.375d0*H33

               endif

               hsep(1,ll+1,ispin)=HH11
               hsep(2,ll+1,ispin)=HH12
               hsep(3,ll+1,ispin)=HH22
               hsep(4,ll+1,ispin)=HH13
               hsep(5,ll+1,ispin)=HH23
               hsep(6,ll+1,ispin)=HH33

            enddo
         enddo
 
         write(6,'(i4,t60,a)') lpx ,
     :        'lpx, (Projectors for l=0..lpx)'
         do l=0,lpx
            write(6,*) 'l=',l
            write(6,'(f7.3,t8,6e11.3,t76,a)') r_l(l+1),
     :           (ohsep(i,l+1,1)-hsep(i,l+1,1),i=1,6),
     :           'r_l(),hsep(), '//is(1)
            if (l.gt.0 .and. nspin.eq.2) 
     :           write(6,'(t8,6e11.3,t76,a)') 
     :           (ohsep(i,l+1,2)-hsep(i,l+1,2),i=1,6),
     :           '       hsep(), '//is(2)
         enddo

      endif

c     restote original  hsep()-values
      do ll=0,lpx
         do i=1,min(2*ll+1,nspin)
            do j=1,6
               hsep(j,ll+1,i)= ohsep(j,ll+1,i)
            enddo
         enddo
      enddo


      end


