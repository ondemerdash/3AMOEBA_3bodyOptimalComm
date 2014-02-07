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
      subroutine cutoffs
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
c         vdwcut = 9.0d0
c         chgcut = 9.0d0
c         dplcut = 9.0d0
c         mpolecut = 9.0d0
         vdwcut = big
         chgcut = big
         dplcut = big
         mpolecut = big
      else
         vdwcut = big
         chgcut = big
         dplcut = big
         mpolecut = big
      end if
      ewaldcut = 7.0d0
c
c     set defaults for tapering, Hessian cutoff and neighbor buffer
c
      vdwtaper = 0.90d0
      chgtaper = 0.65d0
      dpltaper = 0.75d0
      mpoletaper = 0.65d0
      hesscut = 0.0d0
      lbuffer = 2.0d0
c
c     set defaults for Ewald sum, tapering style and neighbor method
c
      use_ewald = .false.
      truncate = .false.
      use_lights = .false.
      use_list = .false.
      use_vlist = .false.
      use_clist = .false.
      use_mlist = .false.
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
c         else if (keyword(1:15) .eq. 'EWALD-CUTOFF3B ') then
c            read (string,*,err=10,end=10)  ewaldcut3b
c
c     get values for the tapering style and neighbor method
c
         else if (keyword(1:9) .eq. 'TRUNCATE ') then
            truncate = .true.
         else if (keyword(1:7) .eq. 'LIGHTS ') then
            use_lights = .true.
         else if (keyword(1:14) .eq. 'NEIGHBOR-LIST ') then
            use_list = .true.
            use_vlist = .true.
            use_clist = .true.
            use_mlist = .true.
         else if (keyword(1:9) .eq. 'VDW-LIST ') then
            use_list = .true.
            use_vlist = .true.
         else if (keyword(1:9) .eq. 'CHG-LIST ') then
            use_list = .true.
            use_clist = .true.
         else if (keyword(1:11) .eq. 'MPOLE-LIST ') then
            use_list = .true.
            use_mlist = .true.
c
c     get cutoff for the magnitude of Hessian elements
c
         else if (keyword(1:12) .eq. 'HESS-CUTOFF ') then
            read (string,*,err=10,end=10)  hesscut
c
c     get the cutoff radii for potential energy functions
c
         else if (keyword(1:7) .eq. 'CUTOFF ') then
            read (string,*,err=10,end=10)  value
            vdwcut = value
            chgcut = value
            dplcut = value
            mpolecut = value
            ewaldcut = value
         else if (keyword(1:11) .eq. 'VDW-CUTOFF ') then
            read (string,*,err=10,end=10)  vdwcut
         else if (keyword(1:11) .eq. 'CHG-CUTOFF ') then
            read (string,*,err=10,end=10)  chgcut
         else if (keyword(1:11) .eq. 'DPL-CUTOFF ') then
            read (string,*,err=10,end=10)  dplcut
         else if (keyword(1:13) .eq. 'MPOLE-CUTOFF ') then
            read (string,*,err=10,end=10)  mpolecut
c
c     get distance for initialization of energy switching
c
         else if (keyword(1:6) .eq. 'TAPER ') then
            read (string,*,err=10,end=10)  value
            vdwtaper = value
            chgtaper = value
            dpltaper = value
            mpoletaper = value
         else if (keyword(1:10) .eq. 'VDW-TAPER ') then
            read (string,*,err=10,end=10)  vdwtaper
         else if (keyword(1:10) .eq. 'CHG-TAPER ') then
            read (string,*,err=10,end=10)  chgtaper
         else if (keyword(1:10) .eq. 'DPL-TAPER ') then
            read (string,*,err=10,end=10)  dpltaper
         else if (keyword(1:12) .eq. 'MPOLE-TAPER ') then
            read (string,*,err=10,end=10)  mpoletaper
c
c     get buffer width for pairwise nonbonded neighbor lists
c
         else if (keyword(1:12) .eq. 'LIST-BUFFER ') then
            read (string,*,err=10,end=10)  lbuffer
         end if
   10    continue
      end do
c
c     apply any Ewald cutoff to charge and multipole terms
c
      if (use_ewald) then
c         chgcut = ewaldcut
c         mpolecut = ewaldcut
      end if
c
c     convert any tapering percentages to absolute distances
c
      if (vdwtaper .lt. 1.0d0)  vdwtaper = vdwtaper * vdwcut
      if (chgtaper .lt. 1.0d0)  chgtaper = chgtaper * chgcut
      if (dpltaper .lt. 1.0d0)  dpltaper = dpltaper * dplcut
      if (mpoletaper .lt. 1.0d0)  mpoletaper = mpoletaper * mpolecut
c
c     apply truncation cutoffs if they were requested
c
      if (truncate) then
         vdwtaper = big
         chgtaper = big
         dpltaper = big
         mpoletaper = big
      end if
c
c     set buffer region limits for pairwise neighbor lists
c
      if (use_list) then
         lbuf2 = (0.5d0*lbuffer)**2
         vbuf2 = (vdwcut+lbuffer)**2
         cbuf2 = (chgcut+lbuffer)**2
         mbuf2 = (mpolecut+lbuffer)**2
         vbufx = (vdwcut+2.0d0*lbuffer)**2
         cbufx = (chgcut+2.0d0*lbuffer)**2
         mbufx = (mpolecut+2.0d0*lbuffer)**2
      end if
c
c     perform dynamic allocation of some pointer arrays
c
      if (use_vlist) then
         if (associated(nvlst))  deallocate (nvlst)
         if (associated(vlst))  deallocate (vlst)
         if (associated(xvold))  deallocate (xvold)
         if (associated(yvold))  deallocate (yvold)
         if (associated(zvold))  deallocate (zvold)
         allocate (nvlst(n))
         allocate (vlst(maxvlst,n))
         allocate (xvold(n))
         allocate (yvold(n))
         allocate (zvold(n))
      end if
      if (use_clist .or. use_mlist) then
         if (associated(nelst))  deallocate (nelst)
         if (associated(elst))  deallocate (elst)
         allocate (nelst(n))
         allocate (elst(maxelst,n))
      end if
      if (use_clist) then
         if (associated(xcold))  deallocate (xcold)
         if (associated(ycold))  deallocate (ycold)
         if (associated(zcold))  deallocate (zcold)
         allocate (xcold(n))
         allocate (ycold(n))
         allocate (zcold(n))
      end if
      if (use_mlist) then
         if (associated(xmold))  deallocate (xmold)
         if (associated(ymold))  deallocate (ymold)
         if (associated(zmold))  deallocate (zmold)
         allocate (xmold(n))
         allocate (ymold(n))
         allocate (zmold(n))
      end if
      return
      end

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
      subroutine cutoffs2_serial
      implicit none
      include 'sizes3b.i'
      include 'sizes.i'
      include 'atoms.i'
      include 'atoms3b.i'
      include 'bound.i'
      include 'cutoff.i'
      include 'hescut.i'
      include 'keys.i'
      include 'neigh3b.i'
      include 'molcul.i'
      include 'molcul3b.i'
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
         if (associated(xmolold))  deallocate (xmolold)
         if (associated(ymolold))  deallocate (ymolold)
         if (associated(zmolold))  deallocate (zmolold)
         if (associated(xmolold3))  deallocate (xmolold3)
         if (associated(ymolold3))  deallocate (ymolold3)
         if (associated(zmolold3))  deallocate (zmolold3)
         if (associated(nmollst))  deallocate (nmollst)
         if (associated(mollst))  deallocate (mollst)
         if (associated(nmollst3mod))  deallocate (nmollst3mod)
         if (associated(mollst3mod))  deallocate (mollst3mod)
c         if (associated(mol3new)) deallocate (mol3new)
         allocate (xmolold(nmol))
         allocate (ymolold(nmol))
         allocate (zmolold(nmol))
         allocate (xmolold3(nmol))
         allocate (ymolold3(nmol))
         allocate (zmolold3(nmol))
         allocate (nmollst(nmol))
         allocate (mollst(100,nmol))
         allocate (nmollst3mod(nmol))
         allocate (mollst3mod(800,nmol))
c         allocate (mol3new(nmol))
      return
      end

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
      subroutine cutoffs2_parallel
      implicit none
      include 'sizes3b.i'
      include 'sizes.i'
      include 'atoms.i'
      include 'atoms3b.i'
      include 'bound.i'
      include 'cutoff.i'
      include 'hescut.i'
      include 'keys.i'
      include 'neigh3b.i'
      include 'molcul.i'
      include 'molcul3b.i'
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

c         if (associated(xmolold))  deallocate (xmolold)
c         if (associated(ymolold))  deallocate (ymolold)
c         if (associated(zmolold))  deallocate (zmolold)
c         if (associated(xmolold3))  deallocate (xmolold3)
c         if (associated(ymolold3))  deallocate (ymolold3)
c         if (associated(zmolold3))  deallocate (zmolold3)
         if (associated(nmollst))  deallocate (nmollst)
         if (associated(mollst))  deallocate (mollst)
c         if (associated(nmollst3mod))  deallocate (nmollst3mod)
c         if (associated(mollst3mod))  deallocate (mollst3mod)
c         if (associated(mol3new)) deallocate (mol3new)
c         allocate (xmolold(nmol))
c         allocate (ymolold(nmol))
c         allocate (zmolold(nmol))
c         allocate (xmolold3(nmol))
c         allocate (ymolold3(nmol))
c         allocate (zmolold3(nmol))
         allocate (nmollst(nmol))
         allocate (mollst(100,nmol))
c         allocate (nmollst3mod(nmol))
c         allocate (mollst3mod(800,nmol))
c         allocate (mol3new(nmol))
      return
      end

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
      subroutine cutoffs_serial
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
         vdwcut = 9.0d0
         chgcut = 9.0d0
         dplcut = 9.0d0
c         mpolecut = 9.0d0
c         vdwcut = big
c         chgcut = big
c         dplcut = big
         mpolecut = big
      else
         vdwcut = big
         chgcut = big
         dplcut = big
         mpolecut = big
      end if
      ewaldcut = 7.0d0
c
c     set defaults for tapering, Hessian cutoff and neighbor buffer
c
      vdwtaper = 0.90d0
      chgtaper = 0.65d0
      dpltaper = 0.75d0
      mpoletaper = 0.65d0
      hesscut = 0.0d0
      lbuffer = 2.0d0
c
c     set defaults for Ewald sum, tapering style and neighbor method
c
      use_ewald = .false.
      truncate = .false.
      use_lights = .false.
      use_list = .false.
      use_vlist = .false.
      use_clist = .false.
      use_mlist = .false.
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
         else if (keyword(1:8) .eq. 'EWALD3B ') then
            use_ewald3b = .true.
         else if (keyword(1:13) .eq. 'EWALD-CUTOFF ') then
            read (string,*,err=10,end=10)  ewaldcut
         else if (keyword(1:15) .eq. 'EWALD-CUTOFF3B ') then
            read (string,*,err=10,end=10)  ewaldcut3b
c
c     get values for the tapering style and neighbor method
c
         else if (keyword(1:9) .eq. 'TRUNCATE ') then
            truncate = .true.
         else if (keyword(1:7) .eq. 'LIGHTS ') then
            use_lights = .true.
         else if (keyword(1:14) .eq. 'NEIGHBOR-LIST ') then
            use_list = .true.
            use_vlist = .true.
            use_clist = .true.
            use_mlist = .true.
         else if (keyword(1:9) .eq. 'VDW-LIST ') then
            use_list = .true.
            use_vlist = .true.
         else if (keyword(1:9) .eq. 'CHG-LIST ') then
            use_list = .true.
            use_clist = .true.
         else if (keyword(1:11) .eq. 'MPOLE-LIST ') then
            use_list = .true.
            use_mlist = .true.
c
c     get cutoff for the magnitude of Hessian elements
c
         else if (keyword(1:12) .eq. 'HESS-CUTOFF ') then
            read (string,*,err=10,end=10)  hesscut
c
c     get the cutoff radii for potential energy functions
c
         else if (keyword(1:7) .eq. 'CUTOFF ') then
            read (string,*,err=10,end=10)  value
            vdwcut = value
            chgcut = value
            dplcut = value
            mpolecut = value
            ewaldcut = value
         else if (keyword(1:11) .eq. 'VDW-CUTOFF ') then
            read (string,*,err=10,end=10)  vdwcut
         else if (keyword(1:11) .eq. 'CHG-CUTOFF ') then
            read (string,*,err=10,end=10)  chgcut
         else if (keyword(1:11) .eq. 'DPL-CUTOFF ') then
            read (string,*,err=10,end=10)  dplcut
         else if (keyword(1:13) .eq. 'MPOLE-CUTOFF ') then
            read (string,*,err=10,end=10)  mpolecut
c
c     get distance for initialization of energy switching
c
         else if (keyword(1:6) .eq. 'TAPER ') then
            read (string,*,err=10,end=10)  value
            vdwtaper = value
            chgtaper = value
            dpltaper = value
            mpoletaper = value
         else if (keyword(1:10) .eq. 'VDW-TAPER ') then
            read (string,*,err=10,end=10)  vdwtaper
         else if (keyword(1:10) .eq. 'CHG-TAPER ') then
            read (string,*,err=10,end=10)  chgtaper
         else if (keyword(1:10) .eq. 'DPL-TAPER ') then
            read (string,*,err=10,end=10)  dpltaper
         else if (keyword(1:12) .eq. 'MPOLE-TAPER ') then
            read (string,*,err=10,end=10)  mpoletaper
c
c     get buffer width for pairwise nonbonded neighbor lists
c
         else if (keyword(1:12) .eq. 'LIST-BUFFER ') then
            read (string,*,err=10,end=10)  lbuffer
         end if
   10    continue
      end do
c
c     apply any Ewald cutoff to charge and multipole terms
c
      if (use_ewald) then
c         chgcut = ewaldcut
c         mpolecut = ewaldcut
      end if
c
c     convert any tapering percentages to absolute distances
c
      if (vdwtaper .lt. 1.0d0)  vdwtaper = vdwtaper * vdwcut
      if (chgtaper .lt. 1.0d0)  chgtaper = chgtaper * chgcut
      if (dpltaper .lt. 1.0d0)  dpltaper = dpltaper * dplcut
      if (mpoletaper .lt. 1.0d0)  mpoletaper = mpoletaper * mpolecut
c
c     apply truncation cutoffs if they were requested
c
      if (truncate) then
         vdwtaper = big
         chgtaper = big
         dpltaper = big
         mpoletaper = big
      end if
c
c     set buffer region limits for pairwise neighbor lists
c
      if (use_list) then
         lbuf2 = (0.5d0*lbuffer)**2
         vbuf2 = (vdwcut+lbuffer)**2
         cbuf2 = (chgcut+lbuffer)**2
         mbuf2 = (mpolecut+lbuffer)**2
         vbufx = (vdwcut+2.0d0*lbuffer)**2
         cbufx = (chgcut+2.0d0*lbuffer)**2
         mbufx = (mpolecut+2.0d0*lbuffer)**2
      end if
c
c     perform dynamic allocation of some pointer arrays
c
      if (use_vlist) then
         if (associated(nvlst))  deallocate (nvlst)
         if (associated(vlst))  deallocate (vlst)
         if (associated(xvold))  deallocate (xvold)
         if (associated(yvold))  deallocate (yvold)
         if (associated(zvold))  deallocate (zvold)
         allocate (nvlst(n))
         allocate (vlst(maxvlst,n))
         allocate (xvold(n))
         allocate (yvold(n))
         allocate (zvold(n))
      end if
      if (use_clist .or. use_mlist) then
         if (associated(nelst))  deallocate (nelst)
         if (associated(elst))  deallocate (elst)
         allocate (nelst(n))
         allocate (elst(maxelst,n))
      end if
      if (use_clist) then
         if (associated(xcold))  deallocate (xcold)
         if (associated(ycold))  deallocate (ycold)
         if (associated(zcold))  deallocate (zcold)
         allocate (xcold(n))
         allocate (ycold(n))
         allocate (zcold(n))
      end if
      if (use_mlist) then
         if (associated(xmold))  deallocate (xmold)
         if (associated(ymold))  deallocate (ymold)
         if (associated(zmold))  deallocate (zmold)
         allocate (xmold(n))
         allocate (ymold(n))
         allocate (zmold(n))
      end if
      return
      end

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
         vdwcut = 9.0d0
      else
         mpolecut = big
      end if
      ewaldcut = 7.0d0
c
c     set defaults for tapering, Hessian cutoff and neighbor buffer
c
      vdwtaper = 0.90d0
      chgtaper = 0.65d0
      dpltaper = 0.75d0
      mpoletaper = 0.65d0
      hesscut = 0.0d0
      lbuffer = 2.0d0

c
c     set defaults for Ewald sum, tapering style and neighbor method
c
      use_ewald = .false.
      truncate = .false.
      use_lights = .false.
      use_list = .false.
      use_vlist = .false.
      use_clist = .false.
      use_mlist = .false.
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

         else if (keyword(1:9) .eq. 'TRUNCATE ') then
            truncate = .true.
         else if (keyword(1:7) .eq. 'LIGHTS ') then
            use_lights = .true.
         else if (keyword(1:14) .eq. 'NEIGHBOR-LIST ') then
            use_list = .true.
            use_vlist = .true.
            use_clist = .true.
            use_mlist = .true.
         else if (keyword(1:9) .eq. 'VDW-LIST ') then
            use_list = .true.
            use_vlist = .true.
         else if (keyword(1:9) .eq. 'CHG-LIST ') then
            use_list = .true.
            use_clist = .true.
         else if (keyword(1:11) .eq. 'MPOLE-LIST ') then
            use_list = .true.
            use_mlist = .true.

c
c     get cutoff for the magnitude of Hessian elements
c

         else if (keyword(1:12) .eq. 'HESS-CUTOFF ') then
            read (string,*,err=10,end=10)  hesscut

c
c     get the cutoff radii for potential energy functions
c

         else if (keyword(1:7) .eq. 'CUTOFF ') then
            read (string,*,err=10,end=10)  value
            vdwcut = value
            chgcut = value
            dplcut = value
            mpolecut = value
            ewaldcut = value
         else if (keyword(1:11) .eq. 'VDW-CUTOFF ') then
            read (string,*,err=10,end=10)  vdwcut
         else if (keyword(1:11) .eq. 'CHG-CUTOFF ') then
            read (string,*,err=10,end=10)  chgcut
         else if (keyword(1:11) .eq. 'DPL-CUTOFF ') then
            read (string,*,err=10,end=10)  dplcut
         else if (keyword(1:13) .eq. 'MPOLE-CUTOFF ') then
            read (string,*,err=10,end=10)  mpolecut

c
c     get distance for initialization of energy switching
c

         else if (keyword(1:6) .eq. 'TAPER ') then
            read (string,*,err=10,end=10)  value
            vdwtaper = value
            chgtaper = value
            dpltaper = value
            mpoletaper = value
         else if (keyword(1:10) .eq. 'VDW-TAPER ') then
            read (string,*,err=10,end=10)  vdwtaper
         else if (keyword(1:10) .eq. 'CHG-TAPER ') then
            read (string,*,err=10,end=10)  chgtaper
         else if (keyword(1:10) .eq. 'DPL-TAPER ') then
            read (string,*,err=10,end=10)  dpltaper
         else if (keyword(1:12) .eq. 'MPOLE-TAPER ') then
            read (string,*,err=10,end=10)  mpoletaper

c
c     get buffer width for pairwise nonbonded neighbor lists
c

         else if (keyword(1:12) .eq. 'LIST-BUFFER ') then
            read (string,*,err=10,end=10)  lbuffer
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

      if (vdwtaper .lt. 1.0d0)  vdwtaper = vdwtaper * vdwcut
      if (chgtaper .lt. 1.0d0)  chgtaper = chgtaper * chgcut
      if (dpltaper .lt. 1.0d0)  dpltaper = dpltaper * dplcut
      if (mpoletaper .lt. 1.0d0)  mpoletaper = mpoletaper * mpolecut

c
c     apply truncation cutoffs if they were requested
c

      if (truncate) then
         vdwtaper = big
         chgtaper = big
         dpltaper = big
         mpoletaper = big
      end if

c
c     set buffer region limits for pairwise neighbor lists
c

      if (use_list) then
         lbuf2 = (0.5d0*lbuffer)**2
         vbuf2 = (vdwcut+lbuffer)**2
         cbuf2 = (chgcut+lbuffer)**2
         mbuf2 = (mpolecut+lbuffer)**2
         vbufx = (vdwcut+2.0d0*lbuffer)**2
         cbufx = (chgcut+2.0d0*lbuffer)**2
         mbufx = (mpolecut+2.0d0*lbuffer)**2
      end if

c
c     perform dynamic allocation of some pointer arrays
c

      if (use_vlist) then
         if (associated(nvlst))  deallocate (nvlst)
         if (associated(vlst))  deallocate (vlst)
         if (associated(xvold))  deallocate (xvold)
         if (associated(yvold))  deallocate (yvold)
         if (associated(zvold))  deallocate (zvold)
         allocate (nvlst(n))
         allocate (vlst(maxvlst,n))
         allocate (xvold(n))
         allocate (yvold(n))
         allocate (zvold(n))
      end if
      if (use_clist .or. use_mlist) then
         if (associated(nelst))  deallocate (nelst)
         if (associated(elst))  deallocate (elst)
         allocate (nelst(n))
         allocate (elst(maxelst,n))
      end if
      if (use_clist) then
         if (associated(xcold))  deallocate (xcold)
         if (associated(ycold))  deallocate (ycold)
         if (associated(zcold))  deallocate (zcold)
         allocate (xcold(n))
         allocate (ycold(n))
         allocate (zcold(n))
      end if
      if (use_mlist) then
         if (associated(xmold))  deallocate (xmold)
         if (associated(ymold))  deallocate (ymold)
         if (associated(zmold))  deallocate (zmold)
         allocate (xmold(n))
         allocate (ymold(n))
         allocate (zmold(n))
      end if
      return
      end

