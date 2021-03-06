c
c
c     #############################################################
c     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  routines below find the respective numbers of atoms in    ##
c     ##  for the interactions between combinations of molecules    ##
c     ##                                                            ##
c     ################################################################
c
c     #############################################################
c     ##                                                         ##
c     ##       subroutine combo1  --  combinations of atoms      ##
c     ##                                                         ##
c     #############################################################
c
c
c     "combo1" finds the index of atoms for which the interactions
c     are being evaluated in empole3
c
c
      subroutine combo1
      implicit none
      include 'sizes3b.i'
      include 'sizes.i'
      include 'mpole3b.i'
      include 'mpole.i'
      include 'mpole3b_3.i'
      include 'mpole3b_2.i'
      include 'molcul.i'
      include 'molcul3b.i'
      include 'combo.i'
      integer i,index
      integer hatom1,latom1
      integer hpole1,lpole1

c
c Find the pole numbers and atom indices
c for each of the 1 body interaction
c
      latom1 = imol(1,moli1)
      hatom1 = imol(2,moli1)
      do i = 1, npole
         if (ipole(i) .eq. latom1) lpole1 = i
         if (ipole(i) .eq. hatom1) hpole1 = i
      end do
      np1 = (hpole1 - lpole1) + 1
      na1 = (hatom1 - latom1) + 1
      index = 0
      do i = 1, np1
         pnum(i) = lpole1 + index
         index = index + 1
      end do
      index = 0
      do i = 1, na1
         anum(i) = latom1 + index
         index = index + 1
      end do
      npole3b = np1
      numatom(moli) = na1
      return
      end

c     #############################################################
c     ##                                                         ##
c     ##       subroutine combo2  --  combinations of atoms      ##
c     ##                                                         ##
c     #############################################################
c
c
c     "combo2" finds the index of atoms for which the interactions
c     are being evaluated in empole3
c
      subroutine combo2
      implicit none
      include 'sizes3b.i'
      include 'sizes.i'
      include 'mpole3b.i'
      include 'mpole.i'
      include 'mpole3b_3.i'
      include 'mpole3b_2.i'
      include 'molcul.i'
      include 'molcul3b.i'
      include 'combo.i'
      integer i,index
      integer hatom1,hatom2
      integer hpole1,hpole2
      integer latom1,latom2
      integer lpole1,lpole2


      latom2 = imol(1,moli2)
      hatom2 = imol(2,moli2)
      do i = 1, npole
          if (ipole(i) .eq. latom2) lpole2 = i
          if (ipole(i) .eq. hatom2) hpole2 = i
      end do
      np2 = np1 + (hpole2 - lpole2) + 1
      na2 = na1 + (hatom2 - latom2) + 1
      index = 0
      do i = (np1+1), np2
          pnum(i) = lpole2 + index
          index = index + 1
      end do
      index = 0
      do i = (na1+1), na2
          anum(i) = latom2 + index
          index = index + 1
      end do
      index = 0
      do i = (na1+1), na2
         anum(i) = imol(1,moli2)+index
         index = index + 1
      end do
      npole3b = np2
      numatom(moli) = na2
      return
      end



c     #############################################################
c     ##                                                         ##
c     ##       subroutine combo3  --  combinations of atoms      ##
c     ##                                                         ##
c     #############################################################
c
c
c     "combo3" finds the index of atoms for which the interactions
c     are being evaluated in empole3
c
      subroutine combo3
      implicit none
      include 'sizes3b.i'
      include 'sizes.i'
      include 'mpole3b.i'
      include 'mpole.i'
      include 'mpole3b_3.i'
      include 'mpole3b_2.i'
      include 'molcul.i'
      include 'molcul3b.i'
      include 'combo.i'
      integer i,index
      integer hatom3, latom3
      integer hpole3, lpole3

  
c
c Get atomic indices of molecule 3
c
      latom3 = imol(1,moli3)
      hatom3 = imol(2,moli3)
      do i = 1, npole
          if (ipole(i) .eq. latom3) lpole3 = i
          if (ipole(i) .eq. hatom3) hpole3 = i
      end do
      np3 = np2 + (hpole3 - lpole3) + 1
      na3 = na2 + (hatom3 - latom3) + 1
      index = 0
      do i = (np2+1), np3
          pnum(i) = lpole3 + index
          index = index + 1
      end do
      index = 0
      do i = (na2+1), na3
          anum(i) = latom3 + index
          index = index + 1
      end do
      index = 0
      do i = (na2+1), na3
         anum(i) = imol(1,moli3)+index
         index = index + 1
      end do
      npole3b = np3
      numatom(moli) = na3
      return
      end


