c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  ewald.i  --  parameters and options for Ewald summation  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     aewald     Ewald convergence coefficient value (Ang-1)
c     boundary   Ewald boundary condition; none, tinfoil or vacuum
c
c
      real*8 aewaldPerm,aewald3b
      integer jmaxPerm,kmaxPerm,lmaxPerm
      integer jmax3b,kmax3b,lmax3b
      character*7 boundary
      common /ewald/ aewaldPerm,aewald3b,
     &               jmaxPerm,kmaxPerm,
     &               lmaxPerm,jmax3b,kmax3b,lmax3b,boundary
