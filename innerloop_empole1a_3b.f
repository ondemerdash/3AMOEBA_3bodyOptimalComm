c
c     ###############################################################
c     ##                                                           ##
c     ##               Subroutine empole1c_3b                      ##
c     ##  Liam O-Suilleabhain, Omar Demerdash, Teresa Head-Gordon  ##
c     ##                 Spring 2013                               ##
c     ###############################################################
c
c
c     "empole1c_3b" calculates the atomic multipole and dipole
c     polarizability interaction energy using the 3-body approximation
c
c
      subroutine Innerloop1(moli1,ep3bt,virep3bt,dep3bt) 
      implicit none
      include 'sizes3b.i'
      include 'sizes.i'
      include 'energi.i'
      include 'molcul.i'
      include 'molcul3b.i'
      include 'deriv.i'
      include 'mpole3b.i'
      include 'mpole.i'
      include 'mpole3b_3.i'
      include 'mpole3b_2.i'
      include 'virial.i'
      include 'mpif.h'
      include 'neigh3b.i'
      real*8  ep1moli1,ep1moli2,ep2moli12,ep1moli3
      real*8  dep1moli1(3,30),dep1moli2(3,30)
      real*8  dep2moli12(3,30),dep1moli3(3,30)
      integer i,ii,j,l1,i1,i2,i3,k
      integer k1,k2,k5
      real*8 eptemp,deptemp(3,npole)
      real*8 vir1moli1(3,3),vir2moli12(3,3),vir1moli3(3,3)
      real*8 vir1moli2(3,3),virtemp(3,3)
      integer pnum(30),npole3b,moli1,moli2,moli3
      real*8 ep3bt,dep3bt(3,*),virep3bt(3,3)
      logical do2
                ep3bt=0.0d0
c                ep3b=0.0d0
                do i = 1, npole
                   do j = 1, 3
                      dep3bt(j,i) = 0.0d0
c                      dep3b(j,i) = 0.0d0
                   end do
                end do

                do i=1,3
                   do j=1,3
c                      virep3b(i,j)=0.0d0
                      virep3bt(i,j)=0.0d0
                   end do
                end do

      do i = 1, 30
        do j = 1, 3
          dep1moli1(j,i)=0.0d0
          dep2moli12(j,i)=0.0d0
          dep1moli2(j,i)=0.0d0
        end do
      end do

      do i=1,3
         do j=1,3
           vir1moli1(i,j)=0.0d0
           vir2moli12(i,j)=0.0d0
           vir1moli2(i,j)=0.0d0
         end do
      end do 

        do k1=1,nmollst(moli1)
          moli2=mollst(k1,moli1)
          
          pnum(1)=imol(1,moli1)
          pnum(2)=imol(1,moli1)+1
          pnum(3)=imol(2,moli1)
          npole3b=6
          pnum(4)=imol(1,moli2)
          pnum(5)=imol(1,moli2)+1
          pnum(6)=imol(2,moli2)
        call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)

          ep3bt = ep3bt + eptemp
          ep2moli12=eptemp
          do l1 = 1, npole3b
            i = pnum(l1)
            do j = 1, 3
              dep3bt(j,i) = dep3bt(j,i)+deptemp(j,i)
              dep2moli12(j,l1)=deptemp(j,i)
            end do
          end do
          do i=1,3
             do j=1,3
                virep3bt(i,j)=virep3bt(i,j)+virtemp(i,j)
                vir2moli12(i,j)=virtemp(i,j)
             end do
          end do

c          do k2=1,nmollst3mod2(moli2,moli1) 
c            done(moli2)=.true.
          do k2=1,nmollst(moli2)
             moli3=mollst(k2,moli2)
c            moli3=mollst3mod2(k2,moli2,moli1)

            npole3b=9
            pnum(1)=imol(1,moli1)
            pnum(2)=imol(1,moli1)+1
            pnum(3)=imol(2,moli1)
            pnum(4)=imol(1,moli2)
            pnum(5)=imol(1,moli2)+1
            pnum(6)=imol(2,moli2)
            pnum(7)=imol(1,moli3)
            pnum(8)=imol(1,moli3)+1
            pnum(9)=imol(2,moli3)
            call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)
            ep3bt = ep3bt + eptemp
            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
                dep3bt(j,i) = dep3bt(j,i)+deptemp(j,i)
              end do
            end do
            do i=1,3
               do j=1,3
                 virep3bt(i,j)=virep3bt(i,j)+virtemp(i,j) 
               end do
            end do
            
            npole3b=6
            ep3bt = ep3bt - ep2moli12
            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
                dep3bt(j,i) = dep3bt(j,i)-dep2moli12(j,l1)
              end do
            end do
            do i=1,3
               do j=1,3
                virep3bt(i,j)=virep3bt(i,j)-vir2moli12(i,j)
               end do
            end do

            npole3b=6
            pnum(1)=imol(1,moli2)
            pnum(2)=imol(1,moli2)+1
            pnum(3)=imol(2,moli2)
            pnum(4)=imol(1,moli3)
            pnum(5)=imol(1,moli3)+1
            pnum(6)=imol(2,moli3)
        call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)
            ep3bt = ep3bt - eptemp
            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
                dep3bt(j,i) = dep3bt(j,i)-deptemp(j,i)
              end do
            end do
            do i=1,3
               do j=1,3
                virep3bt(i,j)=virep3bt(i,j)-virtemp(i,j)
               end do
            end do

            do2=.false.
            do k5=1,nmollst(moli1)
               if(mollst(k5,moli1).eq.moli3) then
                 do2=.true.
                 goto 31
               end if
            end do

   31 continue
          if(do2.eq..true.) then
            npole3b=6
            pnum(1)=imol(1,moli1)
            pnum(2)=imol(1,moli1)+1
            pnum(3)=imol(2,moli1)
            pnum(4)=imol(1,moli3)
            pnum(5)=imol(1,moli3)+1
            pnum(6)=imol(2,moli3)
            call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)
            ep3bt = ep3bt - eptemp
            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
                dep3bt(j,i) = dep3bt(j,i)-deptemp(j,i)
              end do
            end do
            do i=1,3
               do j=1,3
                virep3bt(i,j)=virep3bt(i,j)-virtemp(i,j)
               end do
            end do
          end if


          end do
        end do 
 
c      print*,"InnerloopMoli1= ep3bt=",moli1,ep3bt
      return
      end


      subroutine getlocalsum(moli1,localsum)
      integer localsum,moli1

      localsum=0
      localsum=localsum+1
      print*,"Moli1 localsum",moli1,localsum
      return
      end

c     #############################################################
c     ##                                                         ##
c     ##  subroutine Innerloop2                                  ##
c     ##                                                         ##
c     #############################################################
c
c     Innerloop2 calculates the 2-body polarization energy,gradient, and
c     virial for a single iteration using a neighbor list.


      subroutine Innerloop2(moli1,ep3bt,virep3bt,dep3bt)
      implicit none
      include 'sizes3b.i'
c      include 'energi.i'
      include 'molcul3b.i'
c      include 'deriv.i'
c      include 'mpole.i'
      include 'mpole3b_3.i'
c      include 'virial.i'
c      include 'mpif.h'
      include 'neigh3b.i'
      integer i,ii,j,l1,i1,i2,i3,k
      integer k1,k2,k5
c      real*8 eptemp,deptemp(3,npole)
      real*8 eptemp,deptemp(3,30)
      real*8 virtemp(3,3)
      integer pnum(30),npole3b,moli1,moli2,moli3
      real*8 ep3bt,dep3bt(3,*),virep3bt(3,3)
      logical do2
                ep3bt=0.0d0
                do i = 1, npole
                   do j = 1, 3
                      dep3bt(j,i) = 0.0d0
                   end do
                end do

                do i=1,3
                   do j=1,3
                      virep3bt(i,j)=0.0d0
                   end do
                end do
        do k1=1,nmollst(moli1)
          moli2=mollst(k1,moli1)

          pnum(1)=imol(1,moli1)
          pnum(2)=imol(1,moli1)+1
          pnum(3)=imol(2,moli1)
          npole3b=6
          pnum(4)=imol(1,moli2)
          pnum(5)=imol(1,moli2)+1
          pnum(6)=imol(2,moli2)
        call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)
          ep3bt = ep3bt + eptemp
          do l1 = 1, npole3b
            i = pnum(l1)
            do j = 1, 3
c              dep3bt(j,i) = dep3bt(j,i)+deptemp(j,i)
              dep3bt(j,i) = dep3bt(j,i)+deptemp(j,l1)
            end do
          end do
          do i=1,3
             do j=1,3
                virep3bt(i,j)=virep3bt(i,j)+virtemp(i,j)
             end do
          end do
        end do   
      return
      end

c     #############################################################
c     ##                                                         ##
c     ##  subroutine Innerloop3                                  ##
c     ##                                                         ##
c     #############################################################
c
c     Innerloop3 calculates the 3-body polarization energy,gradient, and
c     virial for a single iteration using a neighbor list.

      subroutine Innerloop3(iter3,moli1,ep3bt,virep3bt,dep3bt)
      implicit none
      include 'sizes3b.i'
c      include 'energi.i'
      include 'molcul3b.i'
c      include 'deriv.i'
c      include 'mpole.i'
      include 'mpole3b_3.i'
c      include 'virial.i'
c      include 'mpif.h'
      include 'neigh3b.i'
      include 'neigh3b_recv.i'
      integer i,ii,j,l1,i1,i2,i3,k,iter3
      integer k1,k2,k5
c      real*8 eptemp,deptemp(3,npole)
      real*8 eptemp,deptemp(3,30)
      real*8 virtemp(3,3)
      integer pnum(30),npole3b,moli1,moli2,moli3
      real*8 ep3bt,dep3bt(3,*),virep3bt(3,3)
      logical do2
                ep3bt=0.0d0
                do i = 1, npole
                   do j = 1, 3
                      dep3bt(j,i) = 0.0d0
                   end do
                end do

                do i=1,3
                   do j=1,3
                      virep3bt(i,j)=0.0d0
                   end do
                end do
c       triplecount=0
c        do k1=1,nmollst3mod(moli1),2
c           moli2=mollst3mod(k1,moli1)
c           moli3=mollst3mod(k1+1,moli1)
        do k1=1,nmollst3mod_recv(iter3),2
           moli2=mollst3mod_recv(k1,iter3)
           moli3=mollst3mod_recv(k1+1,iter3)

c           triplecount=triplecount+1
            npole3b=9
            pnum(1)=imol(1,moli1)
            pnum(2)=imol(1,moli1)+1
            pnum(3)=imol(2,moli1)
            pnum(4)=imol(1,moli2)
            pnum(5)=imol(1,moli2)+1
            pnum(6)=imol(2,moli2)
            pnum(7)=imol(1,moli3)
            pnum(8)=imol(1,moli3)+1
            pnum(9)=imol(2,moli3)

            call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)
            ep3bt = ep3bt + eptemp
            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
c                dep3bt(j,i) = dep3bt(j,i)+deptemp(j,i)
                dep3bt(j,i) = dep3bt(j,i)+deptemp(j,l1)
              end do
            end do
            do i=1,3
               do j=1,3
                 virep3bt(i,j)=virep3bt(i,j)+virtemp(i,j)
               end do
            end do


            do2=.false.
            do k5=1,nmollst(moli1)
               if(mollst(k5,moli1).eq.moli2) then
                 do2=.true.
                 goto 33
               end if
            end do

   33 continue
          if(do2.eq..true.) then
            npole3b=6
            pnum(1)=imol(1,moli1)
            pnum(2)=imol(1,moli1)+1
            pnum(3)=imol(2,moli1)
            pnum(4)=imol(1,moli2)
            pnum(5)=imol(1,moli2)+1
            pnum(6)=imol(2,moli2)
            call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)

            ep3bt=ep3bt-eptemp
             do l1 = 1, npole3b
               i = pnum(l1)
               do j = 1, 3
c                dep3bt(j,i) = dep3bt(j,i)-deptemp(j,i)
                 dep3bt(j,i) = dep3bt(j,i)-deptemp(j,l1)
               end do
             end do
              do i=1,3
               do j=1,3
                virep3bt(i,j)=virep3bt(i,j)-virtemp(i,j)
               end do
              end do
          end if

            do2=.false.
            do k5=1,nmollst(moli2)
               if(mollst(k5,moli2).eq.moli3) then
                 do2=.true.
                 goto 34
               end if
            end do

   34 continue
          if(do2.eq..true.) then
            npole3b=6
            pnum(1)=imol(1,moli2)
            pnum(2)=imol(1,moli2)+1
            pnum(3)=imol(2,moli2)
            pnum(4)=imol(1,moli3)
            pnum(5)=imol(1,moli3)+1
            pnum(6)=imol(2,moli3)
            call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)
            ep3bt = ep3bt - eptemp
            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
c                dep3bt(j,i) = dep3bt(j,i)-deptemp(j,i)
                dep3bt(j,i) = dep3bt(j,i)-deptemp(j,l1)
              end do
            end do
            do i=1,3
               do j=1,3
                virep3bt(i,j)=virep3bt(i,j)-virtemp(i,j)
               end do
            end do
          end if

            do2=.false.
            do k5=1,nmollst(moli1)
               if(mollst(k5,moli1).eq.moli3) then
                 do2=.true.
                 goto 31
               end if
            end do

   31 continue
          if(do2.eq..true.) then
            npole3b=6
            pnum(1)=imol(1,moli1)
            pnum(2)=imol(1,moli1)+1
            pnum(3)=imol(2,moli1)
            pnum(4)=imol(1,moli3)
            pnum(5)=imol(1,moli3)+1
            pnum(6)=imol(2,moli3)
           call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)
            ep3bt = ep3bt - eptemp
            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
c                dep3bt(j,i) = dep3bt(j,i)-deptemp(j,i)
                dep3bt(j,i) = dep3bt(j,i)-deptemp(j,l1)
              end do
            end do
            do i=1,3
              do j=1,3
                virep3bt(i,j)=virep3bt(i,j)-virtemp(i,j)
              end do
            end do
          end if

        end do
      return
      end
c     ##  subroutine Innerloop3                                  ##
c     ##                                                         ##
c     #############################################################
c
c     Innerloop3 calculates the 3-body polarization energy,gradient, and
c     virial for a single iteration using a neighbor list.

      subroutine Innerloop3_small(iter3,moli1,ep3bt,virep3bt,dep3bt)
      implicit none
      include 'sizes3b.i'
c      include 'energi.i'
      include 'molcul3b.i'
c      include 'deriv.i'
c      include 'mpole.i'
      include 'mpole3b_3.i'
c      include 'virial.i'
c      include 'mpif.h'
      include 'neigh3b.i'
      include 'neigh3b_recv.i'
      integer i,ii,j,l1,i1,i2,i3,k,iter3
      integer k1,k2,k5
c      real*8 eptemp,deptemp(3,npole)
      real*8 eptemp,deptemp(3,30)
      real*8 virtemp(3,3)
      integer pnum(30),npole3b,moli1,moli2,moli3
      real*8 ep3bt,dep3bt(3,*),virep3bt(3,3)
      logical do2
                ep3bt=0.0d0
                do i = 1, npole
                   do j = 1, 3
                      dep3bt(j,i) = 0.0d0
                   end do
                end do

                do i=1,3
                   do j=1,3
                      virep3bt(i,j)=0.0d0
                   end do
                end do
c       triplecount=0
c        do k1=1,nmollst3mod(moli1),2
c           moli2=mollst3mod(k1,moli1)
c           moli3=mollst3mod(k1+1,moli1)
        do k1=1,nmollst3mod_recv_small(iter3),2
           moli2=mollst3mod_recv_small(k1,iter3)
           moli3=mollst3mod_recv_small(k1+1,iter3)

c           triplecount=triplecount+1
            npole3b=9
            pnum(1)=imol(1,moli1)
            pnum(2)=imol(1,moli1)+1
            pnum(3)=imol(2,moli1)
            pnum(4)=imol(1,moli2)
            pnum(5)=imol(1,moli2)+1
            pnum(6)=imol(2,moli2)
            pnum(7)=imol(1,moli3)
            pnum(8)=imol(1,moli3)+1
            pnum(9)=imol(2,moli3)

            call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)
            ep3bt = ep3bt + eptemp
            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
c                dep3bt(j,i) = dep3bt(j,i)+deptemp(j,i)
                dep3bt(j,i) = dep3bt(j,i)+deptemp(j,l1)
              end do
            end do
            do i=1,3
               do j=1,3
                 virep3bt(i,j)=virep3bt(i,j)+virtemp(i,j)
               end do
            end do


            do2=.false.
            do k5=1,nmollst(moli1)
               if(mollst(k5,moli1).eq.moli2) then
                 do2=.true.
                 goto 33
               end if
            end do

   33 continue
          if(do2.eq..true.) then
            npole3b=6
            pnum(1)=imol(1,moli1)
            pnum(2)=imol(1,moli1)+1
            pnum(3)=imol(2,moli1)
            pnum(4)=imol(1,moli2)
            pnum(5)=imol(1,moli2)+1
            pnum(6)=imol(2,moli2)
            call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)

            ep3bt=ep3bt-eptemp
             do l1 = 1, npole3b
               i = pnum(l1)
               do j = 1, 3
c                dep3bt(j,i) = dep3bt(j,i)-deptemp(j,i)
                 dep3bt(j,i) = dep3bt(j,i)-deptemp(j,l1)
               end do
             end do
              do i=1,3
               do j=1,3
                virep3bt(i,j)=virep3bt(i,j)-virtemp(i,j)
               end do
              end do
          end if

            do2=.false.
            do k5=1,nmollst(moli2)
               if(mollst(k5,moli2).eq.moli3) then
                 do2=.true.
                 goto 34
               end if
            end do

   34 continue
          if(do2.eq..true.) then
            npole3b=6
            pnum(1)=imol(1,moli2)
            pnum(2)=imol(1,moli2)+1
            pnum(3)=imol(2,moli2)
            pnum(4)=imol(1,moli3)
            pnum(5)=imol(1,moli3)+1
            pnum(6)=imol(2,moli3)
            call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)
            ep3bt = ep3bt - eptemp
            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
c                dep3bt(j,i) = dep3bt(j,i)-deptemp(j,i)
                dep3bt(j,i) = dep3bt(j,i)-deptemp(j,l1)
              end do
            end do
            do i=1,3
               do j=1,3
                virep3bt(i,j)=virep3bt(i,j)-virtemp(i,j)
               end do
            end do
          end if

            do2=.false.
            do k5=1,nmollst(moli1)
               if(mollst(k5,moli1).eq.moli3) then
                 do2=.true.
                 goto 31
               end if
            end do

   31 continue
          if(do2.eq..true.) then
            npole3b=6
            pnum(1)=imol(1,moli1)
            pnum(2)=imol(1,moli1)+1
            pnum(3)=imol(2,moli1)
            pnum(4)=imol(1,moli3)
            pnum(5)=imol(1,moli3)+1
            pnum(6)=imol(2,moli3)
           call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)
            ep3bt = ep3bt - eptemp
            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
c                dep3bt(j,i) = dep3bt(j,i)-deptemp(j,i)
                dep3bt(j,i) = dep3bt(j,i)-deptemp(j,l1)
              end do
            end do
            do i=1,3
              do j=1,3
                virep3bt(i,j)=virep3bt(i,j)-virtemp(i,j)
              end do
            end do
          end if

        end do
      return
      end
