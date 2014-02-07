c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine cutoffs  --  set distance and Hessian cutoffs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "cutoffs" initializes and stores spherical energy cutoff
c     distance windows, Hessian element and Ewald sum cutoffs,
c     and the pairwise neighbor generation method
c
c
      subroutine cutoffs_parallel
      implicit none
      include 'sizes3b.i'
      include 'sizes.i'
      include 'atoms.i'
      include 'atoms3b.i'
      include 'bound.i'
      include 'cutoff.i'
      include 'hescut.i'
      include 'keys.i'
      include 'neigh.i'
      integer i,next
      real*8 big,value
      logical truncate
      character*20 keyword
      character*120 record
      character*120 string
c
c
c     set defaults for spherical energy cutoff distances
c
      big = 1.0d12
      if (use_bounds) then
c         mpolecut = 9.0d0
         mpolecut = big
      else
         mpolecut = big
      end if
      ewaldcut = 7.0d0
c
c     set defaults for tapering, Hessian cutoff and neighbor buffer
c
c
c     set defaults for Ewald sum, tapering style and neighbor method
c
      use_ewald = .false.

c
c     search the keywords for various cutoff parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
c
c     get values related to use of Ewald summation
c
         if (keyword(1:6) .eq. 'EWALD ') then
            use_ewald = .true.
         else if (keyword(1:13) .eq. 'EWALD-CUTOFF ') then
            read (string,*,err=10,end=10)  ewaldcut
         else if (keyword(1:15) .eq. 'EWALD-CUTOFF3B ') then
            read (string,*,err=10,end=10)  ewaldcut3b
c
c     get values for the tapering style and neighbor method
c

c         else if (keyword(1:9) .eq. 'TRUNCATE ') then
c            truncate = .true.
c         else if (keyword(1:7) .eq. 'LIGHTS ') then
c            use_lights = .true.
c         else if (keyword(1:14) .eq. 'NEIGHBOR-LIST ') then
c            use_list = .true.
c            use_vlist = .true.
c            use_clist = .true.
c            use_mlist = .true.
c         else if (keyword(1:9) .eq. 'VDW-LIST ') then
c            use_list = .true.
c            use_vlist = .true.
c         else if (keyword(1:9) .eq. 'CHG-LIST ') then
c            use_list = .true.
c            use_clist = .true.
c         else if (keyword(1:11) .eq. 'MPOLE-LIST ') then
c            use_list = .true.
c            use_mlist = .true.

c
c     get cutoff for the magnitude of Hessian elements
c

c         else if (keyword(1:12) .eq. 'HESS-CUTOFF ') then
c            read (string,*,err=10,end=10)  hesscut

c
c     get the cutoff radii for potential energy functions
c

c         else if (keyword(1:7) .eq. 'CUTOFF ') then
c            read (string,*,err=10,end=10)  value
c            vdwcut = value
c            chgcut = value
c            dplcut = value
c            mpolecut = value
c            ewaldcut = value
c         else if (keyword(1:11) .eq. 'VDW-CUTOFF ') then
c            read (string,*,err=10,end=10)  vdwcut
c         else if (keyword(1:11) .eq. 'CHG-CUTOFF ') then
c            read (string,*,err=10,end=10)  chgcut
c         else if (keyword(1:11) .eq. 'DPL-CUTOFF ') then
c            read (string,*,err=10,end=10)  dplcut
c         else if (keyword(1:13) .eq. 'MPOLE-CUTOFF ') then
c            read (string,*,err=10,end=10)  mpolecut

c
c     get distance for initialization of energy switching
c

c         else if (keyword(1:6) .eq. 'TAPER ') then
c            read (string,*,err=10,end=10)  value
c            vdwtaper = value
c            chgtaper = value
c            dpltaper = value
c            mpoletaper = value
c         else if (keyword(1:10) .eq. 'VDW-TAPER ') then
c            read (string,*,err=10,end=10)  vdwtaper
c         else if (keyword(1:10) .eq. 'CHG-TAPER ') then
c            read (string,*,err=10,end=10)  chgtaper
c         else if (keyword(1:10) .eq. 'DPL-TAPER ') then
c            read (string,*,err=10,end=10)  dpltaper
c         else if (keyword(1:12) .eq. 'MPOLE-TAPER ') then
c            read (string,*,err=10,end=10)  mpoletaper

c
c     get buffer width for pairwise nonbonded neighbor lists
c

c         else if (keyword(1:12) .eq. 'LIST-BUFFER ') then
c            read (string,*,err=10,end=10)  lbuffer
         end if
   10    continue
      end do
c
c     apply any Ewald cutoff to charge and multipole terms
c
c      if (use_ewald) then
c         chgcut = ewaldcut
c         mpolecut = ewaldcut
c      end if
c
c     convert any tapering percentages to absolute distances
c

c      if (vdwtaper .lt. 1.0d0)  vdwtaper = vdwtaper * vdwcut
c      if (chgtaper .lt. 1.0d0)  chgtaper = chgtaper * chgcut
c      if (dpltaper .lt. 1.0d0)  dpltaper = dpltaper * dplcut
c      if (mpoletaper .lt. 1.0d0)  mpoletaper = mpoletaper * mpolecut

c
c     apply truncation cutoffs if they were requested
c

c      if (truncate) then
c         vdwtaper = big
c         chgtaper = big
c         dpltaper = big
c         mpoletaper = big
c      end if

c
c     set buffer region limits for pairwise neighbor lists
c

c      if (use_list) then
c         lbuf2 = (0.5d0*lbuffer)**2
c         vbuf2 = (vdwcut+lbuffer)**2
c         cbuf2 = (chgcut+lbuffer)**2
c         mbuf2 = (mpolecut+lbuffer)**2
c         vbufx = (vdwcut+2.0d0*lbuffer)**2
c         cbufx = (chgcut+2.0d0*lbuffer)**2
c         mbufx = (mpolecut+2.0d0*lbuffer)**2
c      end if

c
c     perform dynamic allocation of some pointer arrays
c

c      if (use_vlist) then
c         if (associated(nvlst))  deallocate (nvlst)
c         if (associated(vlst))  deallocate (vlst)
c         if (associated(xvold))  deallocate (xvold)
c         if (associated(yvold))  deallocate (yvold)
c         if (associated(zvold))  deallocate (zvold)
c         allocate (nvlst(n))
c         allocate (vlst(maxvlst,n))
c         allocate (xvold(n))
c         allocate (yvold(n))
c         allocate (zvold(n))
c      end if
c      if (use_clist .or. use_mlist) then
c         if (associated(nelst))  deallocate (nelst)
c         if (associated(elst))  deallocate (elst)
c         allocate (nelst(n))
c         allocate (elst(maxelst,n))
c      end if
c      if (use_clist) then
c         if (associated(xcold))  deallocate (xcold)
c         if (associated(ycold))  deallocate (ycold)
c         if (associated(zcold))  deallocate (zcold)
c         allocate (xcold(n))
c         allocate (ycold(n))
c         allocate (zcold(n))
c      end if
c      if (use_mlist) then
c         if (associated(xmold))  deallocate (xmold)
c         if (associated(ymold))  deallocate (ymold)
c         if (associated(zmold))  deallocate (zmold)
c         allocate (xmold(n))
c         allocate (ymold(n))
c         allocate (zmold(n))
c      end if
      return
      end

