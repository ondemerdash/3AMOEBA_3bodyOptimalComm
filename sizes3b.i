c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  sizes.i  --  parameter values to set array dimensions  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "sizes.i" sets values for critical array dimensions used
c     throughout the software; these parameters will fix the size
c     of the largest systems that can be handled; values too large
c     for the computer memory or swap space to accomodate will
c     result in poor performance or outright failure
c
c     parameter:      maximum allowed number of:
c
c     maxatm          atoms in the molecular system
c
c
      integer maxatm,maxval,maxcell
      parameter (maxatm=25000)
      parameter (maxval=8)
      parameter (maxcell=10000)
