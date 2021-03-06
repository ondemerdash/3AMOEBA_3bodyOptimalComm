c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine induce1_3bi  --  double loop multipole analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "induce1_3bi" calculates the direct field from the interaction
c     between atoms in 1-body
c
c
      subroutine induce1_3bi (fdir,field,fieldp)
      implicit none
      include 'boxes.i'
      include 'sizes3b.i'
      include 'sizes.i'
      include 'ewald.i'
      include 'math.i'
      include 'atoms.i'
      include 'atoms3b.i'
      include 'mpole3b.i'
      include 'mpole.i'
      include 'mpole3b_3.i'
      include 'mpole3b_2.i'
      include 'combo.i'
      real*8 term
      real*8 ucell(3)
      real*8 fdir(3,*)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      integer i,ii,j,l1
c
c     get the reciprical space part of the electrostatic field
c
      call udirect1_3b (field)
c
c     get the real space portion of the electrostatic field
c
      do l1 = 1, npole3b
         i = pnum(l1)
         do j = 1, 3
            fieldp(j,i) = field(j,i)
         end do
      end do
      call udirect2a_3b (field,fieldp)
c
c     get the self-energy portion of the electrostatic field
c
      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      do l1 = 1, npole3b
         i = pnum(l1)
         do j = 1, 3
            field(j,i) = field(j,i) + term*rpole(j+1,i)
            fieldp(j,i) = fieldp(j,i) + term*rpole(j+1,i)
         end do
      end do

      if (boundary .eq. 'VACUUM') then
         do i = 1, 3
            ucell(i) = 0.0d0
         end do
         do l1 = 1, npole3b
            i = pnum(l1)
            ii = ipole(i)
            ucell(1) = ucell(1) + rpole(2,i) + rpole(1,i)*x(ii)
            ucell(2) = ucell(2) + rpole(3,i) + rpole(1,i)*y(ii)
            ucell(3) = ucell(3) + rpole(4,i) + rpole(1,i)*z(ii)
         end do
         term = (4.0d0/3.0d0) * pi/volbox
         do l1 = 1, npole3b
            i = pnum(l1)
            do j = 1, 3
               field(j,i) = field(j,i) - term*ucell(j)
               fieldp(j,i) = fieldp(j,i) - term*ucell(j)
            end do
         end do
      end if
      do l1 = 1, npole3b
         i = pnum(l1)
         do j = 1, 3
            fdir(j,i) = field(j,i)
            fd_1b(j,i) = fdir(j,i)
         end do
      end do

      return
      end


c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine induce2_3bi  --  double loop multipole analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "induce2_3bi" calculates the direct field from the interaction
c     between atoms in 1-body
c
c
      subroutine induce2_3bi (fdir,field,fieldp)
      implicit none
      include 'boxes.i'
      include 'sizes3b.i'
      include 'sizes.i'
      include 'ewald.i'
      include 'math.i'
      include 'atoms.i'
      include 'atoms3b.i'
      include 'mpole3b.i'
      include 'mpole.i'
      include 'mpole3b_3.i'
      include 'mpole3b_2.i'
      include 'combo.i'
      real*8 term
      real*8 ucell(3)
      real*8 fdir(3,*)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      integer i,ii,j,l1
c
c     get the reciprical space part of the electrostatic field
c
      call udirect1_3b (field)
c
c     get the real space portion of the electrostatic field
c
      do l1 = 1, npole3b
         i = pnum(l1)
         do j = 1, 3
            fieldp(j,i) = field(j,i)
         end do
      end do
      call udirect2a_3b (field,fieldp)
      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      do l1 = 1, npole3b
         i = pnum(l1)
         do j = 1, 3
            field(j,i) = field(j,i) + term*rpole(j+1,i)
            fieldp(j,i) = fieldp(j,i) + term*rpole(j+1,i)
         end do
      end do

      if (boundary .eq. 'VACUUM') then
         do i = 1, 3
            ucell(i) = 0.0d0
         end do
         do l1 = 1, npole3b
            i = pnum(l1)
            ii = ipole(i)
            ucell(1) = ucell(1) + rpole(2,i) + rpole(1,i)*x(ii)
            ucell(2) = ucell(2) + rpole(3,i) + rpole(1,i)*y(ii)
            ucell(3) = ucell(3) + rpole(4,i) + rpole(1,i)*z(ii)
         end do
         term = (4.0d0/3.0d0) * pi/volbox
         do l1 = 1, npole3b
            i = pnum(l1)
            do j = 1, 3
               field(j,i) = field(j,i) - term*ucell(j)
               fieldp(j,i) = fieldp(j,i) - term*ucell(j)
            end do
         end do
      end if
      do l1 = 1, np1
         i = pnum(l1)
         do j = 1, 3
            fd_2b(j,i,moli2) = field(j,i) - fd_1b(j,i)
         end do
      end do
      do l1 = np1+1, np2
         i = pnum(l1)
         do j = 1, 3
            fd_2b(j,i,moli1) = field(j,i) - fd_1b(j,i)
         end do
      end do
      do l1 = 1, npole3b
         i = pnum(l1)
         do j = 1, 3
            fdir(j,i) = field(j,i)
         end do
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine induce3_3bi  --  double loop multipole analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "induce3_3bi" calculates the direct field from the interaction
c     between atoms in 1-body
c
c
      subroutine induce3_3bi (fdir,field,fieldp)
      implicit none
      include 'boxes.i'
      include 'sizes3b.i'
      include 'sizes.i'
      include 'ewald.i'
      include 'math.i'
      include 'atoms.i'
      include 'atoms3b.i'
      include 'mpole3b.i'
      include 'mpole.i'
      include 'mpole3b_3.i'
      include 'mpole3b_2.i'
      include 'combo.i'
      real*8 term
      real*8 ucell(3)
      real*8 fdir(3,*)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      integer i,ii,j,l1

      do l1 = 1, np1
         i = pnum(l1)
         do j = 1, 3
            fdir(j,i) = fd_2b(j,i,moli2)
     &               + fd_2b(j,i,moli3) + fd_1b(j,i)
         end do
      end do
      do l1 = np1+1, np2
         i = pnum(l1)
         do j = 1, 3
            fdir(j,i) = fd_2b(j,i,moli1)
     &               + fd_2b(j,i,moli3) + fd_1b(j,i)
         end do
      end do
      do l1 = np2+1, np3
         i = pnum(l1)
         do j = 1, 3
            fdir(j,i) = fd_2b(j,i,moli1)
     &               + fd_2b(j,i,moli2) + fd_1b(j,i)
         end do
      end do
      return
      end

