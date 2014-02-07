c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine kewald  --  Ewald sum parameter assignment  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "kewald" assigns particle mesh Ewald parameters and options
c     for a periodic system
c
c
      subroutine kewald3b
      implicit none
      include 'sizes3b.i'
      include 'sizes.i'
      include 'atoms.i'
      include 'atoms3b.i'
      include 'bound.i'
      include 'boxes.i'
      include 'chunks.i'
      include 'cutoff.i'
      include 'ewald.i'
      include 'fft.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'openmp.i'
      include 'pme.i'
      include 'math.i'
      integer maxpower
      parameter (maxpower=54)
      integer i,k,next
      integer ifft1,ifft2,ifft3
      integer multi(maxpower)
      real*8 delta,rmax,dens
      character*20 keyword
      character*120 record
      character*120 string
c
c     PME grid size must be even with factors of only 2, 3 and 5
c
      data multi  /   2,   4,   6,   8,  10,  12,  16,  18,  20,
     &               24,  30,  32,  36,  40,  48,  50,  54,  60,
     &               64,  72,  80,  90,  96, 100, 108, 120, 128,
     &              144, 150, 160, 162, 180, 192, 200, 216, 240,
     &              250, 256, 270, 288, 300, 320, 324, 360, 384,
     &              400, 432, 450, 480, 486, 500, 512, 540, 576 /
c
c
c     default boundary treatment, B-spline order and grid density
c
      ffttyp = 'FFTPACK'
      if (nthread .gt. 1)  ffttyp = 'FFTW'
      boundary = 'TINFOIL'
      bsorder = 5
      dens = 1.2d0
c
c     estimate an optimal value for the Ewald coefficient
c
c      call ewaldcof (aewaldPerm,ewaldcut)
      call ewaldcof (aewald3b,ewaldcut3b)

c      jmaxPerm=int(ewaldcut*aewaldPerm*xbox*aewaldPerm/pi)
c      kmaxPerm=int(ewaldcut*aewaldPerm*ybox*aewaldPerm/pi)
c      lmaxPerm=int(ewaldcut*aewaldPerm*zbox*aewaldPerm/pi)
c      print*,"PermElec vec",jmaxPerm,kmaxPerm,lmaxPerm
c      print*,"PermElec aewald",aewaldPerm

      jmax3b=int(ewaldcut3b*aewald3b*xbox*aewald3b/pi)
      kmax3b=int(ewaldcut3b*aewald3b*ybox*aewald3b/pi)
      lmax3b=int(ewaldcut3b*aewald3b*zbox*aewald3b/pi)
c      print*,"PolElec vec3b",jmax3b,kmax3b,lmax3b
c      print*,"PolElec aewald3b",aewald3b
c
c     set the system extent for nonperiodic Ewald summation
c
      if (.not. use_bounds) then
         call extent (rmax)
         xbox = 2.0d0 * (rmax+ewaldcut)
         ybox = xbox
         zbox = xbox
         alpha = 90.0d0
         beta = 90.0d0
         gamma = 90.0d0
         orthogonal = .true.
         call lattice
         boundary = 'NONE'
         dens = 0.75d0
      end if
c
c     get default grid size from periodic system dimensions
c

c      delta = 1.0d-8
c      ifft1 = int(xbox*dens-delta) + 1
c      ifft2 = int(ybox*dens-delta) + 1
c      ifft3 = int(zbox*dens-delta) + 1

c
c     search keywords for Ewald summation commands
c
      do i = 1, nkey
         record = keyline(i)
         next = 1
         call upcase (record)
         call gettext (record,keyword,next)
         string = record(next:120)
         if (keyword(1:12) .eq. 'FFT-PACKAGE ') then
            call getword (record,ffttyp,next)
         else if (keyword(1:12) .eq. 'EWALD-ALPHA ') then
            read (string,*,err=20)  aewaldPerm
         else if (keyword(1:15) .eq. 'EWALD-BOUNDARY ') then
            boundary = 'VACUUM'
         else if (keyword(1:9) .eq. 'PME-GRID ') then
            ifft1 = 0
            ifft2 = 0
            ifft3 = 0
            read (string,*,err=10,end=10)  ifft1,ifft2,ifft3
   10       continue
            if (ifft2 .eq. 0)  ifft2 = ifft1
            if (ifft3 .eq. 0)  ifft3 = ifft1
         else if (keyword(1:10) .eq. 'PME-ORDER ') then
            read (string,*,err=20)  bsorder
         end if
   20    continue
      end do
c
c     grid size must be even, with prime factors of 2, 3 and 5
c

c      nfft1 = maxfft
c      nfft2 = maxfft
c      nfft3 = maxfft
c      do i = maxpower, 1, -1
c         k = multi(i)
c         if (k .le. maxfft) then
c            if (k .ge. ifft1)  nfft1 = k
c            if (k .ge. ifft2)  nfft2 = k
c            if (k .ge. ifft3)  nfft3 = k
c         end if
c      end do

c
c     set the number of chunks and grid points per chunk
c

c      call getchunk

c
c     check the B-spline order and charge grid dimension
c

c      if (bsorder .gt. maxorder) then
c         write (iout,30)
c   30    format (/,' KEWALD  --  B-Spline Order Too Large;',
c     &              ' Increase MAXORDER')
c         call fatal
c      end if
c      if (max(nfft1,nfft2,nfft3) .gt. maxfft) then
c         write (iout,40)
c   40    format (/,' KEWALD  --  FFT Charge Grid Too Large;',
c     &              ' Increase MAXFFT')
c         call fatal
c      else if (nfft1.lt.ifft1 .or. nfft2.lt.ifft2
c     &             .or. nfft3.lt.ifft3) then
c         write (iout,50)
c   50    format (/,' KEWALD  --  Warning, Small Charge Grid',
c     &              ' may give Poor Accuracy')
c      end if

c
c     perform dynamic allocation of some pointer arrays
c
c      if (associated(thetai1))  deallocate (thetai1)
c      if (associated(thetai2))  deallocate (thetai2)
c      if (associated(thetai3))  deallocate (thetai3)
c      if (associated(qgrid))  deallocate (qgrid)
c      if (associated(qfac))  deallocate (qfac)
c      if (associated(pmetable))  deallocate (pmetable)
c      allocate (thetai1(4,bsorder,n))
c      allocate (thetai2(4,bsorder,n))
c      allocate (thetai3(4,bsorder,n))
c      allocate (qgrid(2,nfft1,nfft2,nfft3))
c      allocate (qfac(nfft1,nfft2,nfft3))
c      allocate (pmetable(n,nchunk))
c
c     initialize the PME arrays that can be precomputed
c

c      call moduli
c      call fftsetup

c
c     print a message listing some of the Ewald parameters
c

c      if (verbose) then
c         write (iout,60)  aewaldPerm,nfft1,nfft2,nfft3,bsorder
c   60    format (/,' Smooth Particle Mesh Ewald Parameters :',
c     &           //,4x,'Ewald Coefficient',6x,'Charge Grid',
c     &              ' Dimensions',6x,'B-Spline Order',
c     &           //,8x,f8.4,11x,3i6,12x,i6)
c      end if
      return
      end
