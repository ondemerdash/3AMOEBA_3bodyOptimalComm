c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine emrecip1  --  mpole Ewald recip energy & derivs  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "emrecip1" evaluates the reciprocal space portion of the particle
c     mesh Ewald summation energy and gradient due to atomic multipole
c     interactions and dipole polarizability
c
c     literature reference:
c
c     C. Sagui, L. G. Pedersen and T. A. Darden, "Towards an Accurate
c     Representation of Electrostatics in Classical Force Fields:
c     Efficient Implementation of Multipolar Interactions in
c     Biomolecular Simulations", Journal of Chemical Physics, 120,
c     73-87 (2004)
c
c     modifications for nonperiodic systems suggested by Tom Darden
c     during May 2007
c
c
      subroutine emrecip1_3b_Perm
      implicit none
c      include 'sizes.i'
c      include 'atoms.i'
c      include 'bound.i'
c      include 'boxes.i'
c      include 'chgpot.i'
c      include 'deriv.i'
c      include 'energi.i'
c      include 'ewald.i'
c      include 'math.i'
c      include 'mpole.i'
c      include 'pme.i'
c      include 'polar.i'
c      include 'polpot.i'
c      include 'potent.i'
c      include 'virial.i'
c      include 'sizes.i'
      include 'sizes3b.i'
      include 'sizes.i'
      include 'atoms.i'
      include 'atoms3b.i'
      include 'bound.i'
      include 'boxes.i'
      include 'chgpot.i'
      include 'chgpot3b.i'
      include 'deriv.i'
      include 'energi.i'
      include 'ewald.i'
      include 'math.i'
c      include 'mpole.i'
      include 'mpole3b.i'
      include 'mpole.i'
      include 'mpole3b_3.i'
      include 'mpole3b_2.i'
      include 'pme.i'
      include 'polar.i'
      include 'polpot.i'
      include 'potent.i'
      include 'virial.i'
      integer i,j,k,ii
      integer j1,j2,j3
      integer k1,k2,k3
      integer m1,m2,m3
      integer ntot,nff
      integer nf1,nf2,nf3
      integer deriv1(10)
      integer deriv2(10)
      integer deriv3(10)
      real*8 e,eterm
      real*8 r1,r2,r3
      real*8 h1,h2,h3
      real*8 f1,f2,f3
      real*8 vxx,vyx,vzx
      real*8 vyy,vzy,vzz
      real*8 volterm,denom
      real*8 hsq,expterm
      real*8 term,pterm
      real*8 vterm,struc2
      real*8 cphim(4),cphid(4)
      real*8 cphip(4)
      real*8 a(3,3),ftc(10,10)
      real*8, allocatable :: frc(:,:)
      real*8, allocatable :: trq(:,:)
      real*8, allocatable :: fuind(:,:)
      real*8, allocatable :: fuinp(:,:)
      real*8, allocatable :: cmp(:,:)
      real*8, allocatable :: fmp(:,:)
      real*8, allocatable :: fphi(:,:)
      real*8, allocatable :: fphid(:,:)
      real*8, allocatable :: fphip(:,:)
      real*8, allocatable :: fphidp(:,:)
      real*8, allocatable :: cphi(:,:)
      real*8, allocatable :: qgrip(:,:,:,:)
      real*8 aewald
c
c     derivative indices into the fphi and fphidp arrays
c
      data deriv1  / 2, 5,  8,  9, 11, 16, 18, 14, 15, 20 /
      data deriv2  / 3, 8,  6, 10, 14, 12, 19, 16, 20, 17 /
      data deriv3  / 4, 9, 10,  7, 15, 17, 13, 20, 18, 19 /
c
c
c     return if the Ewald coefficient is zero
c
      aewald=aewaldPerm
      if (aewald .lt. 1.0d-6)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (frc(3,n))
      allocate (trq(3,npole))
      allocate (fuind(3,npole))
      allocate (fuinp(3,npole))
      allocate (cmp(10,npole))
      allocate (fmp(10,npole))
      allocate (fphi(20,npole))
      allocate (fphid(10,npole))
      allocate (fphip(10,npole))
      allocate (fphidp(20,npole))
      allocate (cphi(10,npole))
c
c     zero out the temporary virial accumulation variables
c
      vxx = 0.0d0
      vyx = 0.0d0
      vzx = 0.0d0
      vyy = 0.0d0
      vzy = 0.0d0
      vzz = 0.0d0
c
c     copy multipole moments and coordinates to local storage
c
      do i = 1, npole
         cmp(1,i) = rpole(1,i)
         cmp(2,i) = rpole(2,i)
         cmp(3,i) = rpole(3,i)
         cmp(4,i) = rpole(4,i)
         cmp(5,i) = rpole(5,i)
         cmp(6,i) = rpole(9,i)
         cmp(7,i) = rpole(13,i)
         cmp(8,i) = 2.0d0 * rpole(6,i)
         cmp(9,i) = 2.0d0 * rpole(7,i)
         cmp(10,i) = 2.0d0 * rpole(10,i)
      end do
c
c     get the fractional to Cartesian transformation matrix
c
      call frac_to_cart (ftc)
c
c     compute the arrays of B-spline coefficients
c
c      if (.not. use_polar) then
         call bspline_fill
         call table_fill
c      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (qgrip(2,nfft1,nfft2,nfft3))
c
c     assign permanent and induced multipoles to PME grid
c     and perform the 3-D FFT forward transformation
c

c      if (use_polar) then
c         do i = 1, npole
c            do j = 2, 4
c               cmp(j,i) = cmp(j,i) + uinp(j-1,i)
c            end do
c         end do
c         call cmp_to_fmp (cmp,fmp)
c         call grid_mpole (fmp)
c         call fftfront
c         do k = 1, nfft3
c            do j = 1, nfft2
c               do i = 1, nfft1
c                  qgrip(1,i,j,k) = qgrid(1,i,j,k)
c                  qgrip(2,i,j,k) = qgrid(2,i,j,k)
c               end do
c            end do
c         end do
c         do i = 1, npole
c            do j = 2, 4
c               cmp(j,i) = cmp(j,i) + uind(j-1,i) - uinp(j-1,i)
c            end do
c         end do
c         call cmp_to_fmp (cmp,fmp)
c         call grid_mpole (fmp)
c         call fftfront
c         do i = 1, npole
c            do j = 2, 4
c               cmp(j,i) = cmp(j,i) - uind(j-1,i)
c            end do
c         end do
c      else
         call cmp_to_fmp (cmp,fmp)
         call grid_mpole (fmp)
         call fftfront
         do k = 1, nfft3
            do j = 1, nfft2
               do i = 1, nfft1
                  qgrip(1,i,j,k) = qgrid(1,i,j,k)
                  qgrip(2,i,j,k) = qgrid(2,i,j,k)
               end do
            end do
         end do
c      end if

c
c     make the scalar summation over reciprocal lattice
c
      ntot = nfft1 * nfft2 * nfft3
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nff = nfft1 * nfft2
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
      do i = 1, ntot-1
         k3 = i/nff + 1
         j = i - (k3-1)*nff
         k2 = j/nfft1 + 1
         k1 = j - (k2-1)*nfft1 + 1
         m1 = k1 - 1
         m2 = k2 - 1
         m3 = k3 - 1
         if (k1 .gt. nf1)  m1 = m1 - nfft1
         if (k2 .gt. nf2)  m2 = m2 - nfft2
         if (k3 .gt. nf3)  m3 = m3 - nfft3
         r1 = dble(m1)
         r2 = dble(m2)
         r3 = dble(m3)
         h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
         h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
         h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
         hsq = h1*h1 + h2*h2 + h3*h3
         term = -pterm * hsq
         expterm = 0.0d0
         if (term .gt. -50.0d0) then
            denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
            expterm = exp(term) / denom
            if (.not. use_bounds) then
               expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
            else if (octahedron) then
               if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
            end if
            struc2 = qgrid(1,k1,k2,k3)*qgrip(1,k1,k2,k3)
     &                  + qgrid(2,k1,k2,k3)*qgrip(2,k1,k2,k3)
            eterm = 0.5d0 * electric * expterm * struc2
            vterm = (2.0d0/hsq) * (1.0d0-term) * eterm
            vxx = vxx + h1*h1*vterm - eterm
            vyx = vyx + h2*h1*vterm
            vzx = vzx + h3*h1*vterm
            vyy = vyy + h2*h2*vterm - eterm
            vzy = vzy + h3*h2*vterm
            vzz = vzz + h3*h3*vterm - eterm
         end if
         qfac(k1,k2,k3) = expterm
      end do
c
c     assign just the induced multipoles to PME grid
c     and perform the 3-D FFT forward transformation
c

c      if (use_polar .and. poltyp.eq.'DIRECT') then
c         do i = 1, npole
c            do j = 1, 10
c               cmp(j,i) = 0.0d0
c            end do
c            do j = 2, 4
c               cmp(j,i) = uinp(j-1,i)
c            end do
c         end do
c         call cmp_to_fmp (cmp,fmp)
c         call grid_mpole (fmp)
c         call fftfront
c         do k = 1, nfft3
c            do j = 1, nfft2
c               do i = 1, nfft1
c                  qgrip(1,i,j,k) = qgrid(1,i,j,k)
c                  qgrip(2,i,j,k) = qgrid(2,i,j,k)
c               end do
c            end do
c         end do
c         do i = 1, npole
c            do j = 2, 4
c               cmp(j,i) = uind(j-1,i)
c            end do
c         end do
c         call cmp_to_fmp (cmp,fmp)
c         call grid_mpole (fmp)
c         call fftfront
c         do i = 1, npole
c            cmp(1,i) = rpole(1,i)
c            cmp(2,i) = rpole(2,i)
c            cmp(3,i) = rpole(3,i)
c            cmp(4,i) = rpole(4,i)
c            cmp(5,i) = rpole(5,i)
c            cmp(6,i) = rpole(9,i)
c            cmp(7,i) = rpole(13,i)
c            cmp(8,i) = 2.0d0 * rpole(6,i)
c            cmp(9,i) = 2.0d0 * rpole(7,i)
c            cmp(10,i) = 2.0d0 * rpole(10,i)
c         end do
c
c     make the scalar summation over reciprocal lattice
c
c         do i = 1, ntot-1
c            k3 = i/nff + 1
c            j = i - (k3-1)*nff
c            k2 = j/nfft1 + 1
c            k1 = j - (k2-1)*nfft1 + 1
c            m1 = k1 - 1
c            m2 = k2 - 1
c            m3 = k3 - 1
c            if (k1 .gt. nf1)  m1 = m1 - nfft1
c            if (k2 .gt. nf2)  m2 = m2 - nfft2
c            if (k3 .gt. nf3)  m3 = m3 - nfft3
c            r1 = dble(m1)
c            r2 = dble(m2)
c            r3 = dble(m3)
c            h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
c            h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
c            h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
c            hsq = h1*h1 + h2*h2 + h3*h3
c            term = -pterm * hsq
c            expterm = 0.0d0
c            if (term .gt. -50.0d0) then
c               denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
c               expterm = exp(term) / denom
c               if (.not. use_bounds) then
c                  expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
c               else if (octahedron) then
c                  if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
c               end if
c               struc2 = qgrid(1,k1,k2,k3)*qgrip(1,k1,k2,k3)
c     &                     + qgrid(2,k1,k2,k3)*qgrip(2,k1,k2,k3)
c               eterm = 0.5d0 * electric * expterm * struc2
c               vterm = (2.0d0/hsq) * (1.0d0-term) * eterm
c               vxx = vxx - h1*h1*vterm + eterm
c               vyx = vyx - h2*h1*vterm
c               vzx = vzx - h3*h1*vterm
c               vyy = vyy - h2*h2*vterm + eterm
c               vzy = vzy - h3*h2*vterm
c               vzz = vzz - h3*h3*vterm + eterm
c            end if
c         end do
c      end if
c
c     perform deallocation of some local arrays
c
      deallocate (qgrip)
c
c     transform permanent multipoles without induced dipoles
c

c      if (use_polar) then
c         call cmp_to_fmp (cmp,fmp)
c         call grid_mpole (fmp)
c         call fftfront
c      end if

c
c     account for the zeroth grid point for a finite system
c
      qfac(1,1,1) = 0.0d0
      if (.not. use_bounds) then
         expterm = 0.5d0 * pi / xbox
         struc2 = qgrid(1,1,1,1)**2 + qgrid(2,1,1,1)**2
         e = 0.5d0 * expterm * struc2
         qfac(1,1,1) = expterm
      end if
c
c     complete the transformation of the PME grid
c
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
               term = qfac(i,j,k)
               qgrid(1,i,j,k) = term * qgrid(1,i,j,k)
               qgrid(2,i,j,k) = term * qgrid(2,i,j,k)
            end do
         end do
      end do
c
c     perform 3-D FFT backward transform and get potential
c
      call fftback
      call fphi_mpole (fphi)
      do i = 1, npole
         do j = 1, 20
            fphi(j,i) = electric * fphi(j,i)
         end do
      end do
      call fphi_to_cphi (fphi,cphi)
c
c     increment the permanent multipole energy and gradient
c
      e = 0.0d0
      do i = 1, npole
         f1 = 0.0d0
         f2 = 0.0d0
         f3 = 0.0d0
         do k = 1, 10
            e = e + fmp(k,i)*fphi(k,i)
            f1 = f1 + fmp(k,i)*fphi(deriv1(k),i)
            f2 = f2 + fmp(k,i)*fphi(deriv2(k),i)
            f3 = f3 + fmp(k,i)*fphi(deriv3(k),i)
         end do
         f1 = dble(nfft1) * f1
         f2 = dble(nfft2) * f2
         f3 = dble(nfft3) * f3
         frc(1,i) = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
         frc(2,i) = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
         frc(3,i) = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
      end do
      e = 0.5d0 * e
      em = em + e
      do i = 1, npole
         ii = ipole(i)
         dem(1,ii) = dem(1,ii) + frc(1,i)
         dem(2,ii) = dem(2,ii) + frc(2,i)
         dem(3,ii) = dem(3,ii) + frc(3,i)
      end do
c
c     distribute torques into the permanent multipole gradient
c
      do i = 1, npole
         trq(1,i) = cmp(4,i)*cphi(3,i) - cmp(3,i)*cphi(4,i)
     &                 + 2.0d0*(cmp(7,i)-cmp(6,i))*cphi(10,i)
     &                 + cmp(9,i)*cphi(8,i) + cmp(10,i)*cphi(6,i)
     &                 - cmp(8,i)*cphi(9,i) - cmp(10,i)*cphi(7,i)
         trq(2,i) = cmp(2,i)*cphi(4,i) - cmp(4,i)*cphi(2,i)
     &                 + 2.0d0*(cmp(5,i)-cmp(7,i))*cphi(9,i)
     &                 + cmp(8,i)*cphi(10,i) + cmp(9,i)*cphi(7,i)
     &                 - cmp(9,i)*cphi(5,i) - cmp(10,i)*cphi(8,i)
         trq(3,i) = cmp(3,i)*cphi(2,i) - cmp(2,i)*cphi(3,i)
     &                 + 2.0d0*(cmp(6,i)-cmp(5,i))*cphi(8,i)
     &                 + cmp(8,i)*cphi(5,i) + cmp(10,i)*cphi(9,i)
     &                 - cmp(8,i)*cphi(6,i) - cmp(9,i)*cphi(10,i)
      end do
      do i = 1, n
         frc(1,i) = 0.0d0
         frc(2,i) = 0.0d0
         frc(3,i) = 0.0d0
      end do
      call torque2 (trq,frc)
      do i = 1, n
         dem(1,i) = dem(1,i) + frc(1,i)
         dem(2,i) = dem(2,i) + frc(2,i)
         dem(3,i) = dem(3,i) + frc(3,i)
      end do
c
c     permanent multipole contribution to the internal virial
c
      do i = 1, npole
         vxx = vxx - cmp(2,i)*cphi(2,i) - 2.0d0*cmp(5,i)*cphi(5,i)
     &             - cmp(8,i)*cphi(8,i) - cmp(9,i)*cphi(9,i)
         vyx = vyx - 0.5d0*(cmp(3,i)*cphi(2,i)+cmp(2,i)*cphi(3,i))
     &             - (cmp(5,i)+cmp(6,i))*cphi(8,i)
     &             - 0.5d0*cmp(8,i)*(cphi(5,i)+cphi(6,i))
     &             - 0.5d0*(cmp(9,i)*cphi(10,i)+cmp(10,i)*cphi(9,i))
         vzx = vzx - 0.5d0*(cmp(4,i)*cphi(2,i)+cmp(2,i)*cphi(4,i))
     &             - (cmp(5,i)+cmp(7,i))*cphi(9,i)
     &             - 0.5d0*cmp(9,i)*(cphi(5,i)+cphi(7,i))
     &             - 0.5d0*(cmp(8,i)*cphi(10,i)+cmp(10,i)*cphi(8,i))
         vyy = vyy - cmp(3,i)*cphi(3,i) - 2.0d0*cmp(6,i)*cphi(6,i)
     &             - cmp(8,i)*cphi(8,i) - cmp(10,i)*cphi(10,i)
         vzy = vzy - 0.5d0*(cmp(4,i)*cphi(3,i)+cmp(3,i)*cphi(4,i))
     &             - (cmp(6,i)+cmp(7,i))*cphi(10,i)
     &             - 0.5d0*cmp(10,i)*(cphi(6,i)+cphi(7,i))
     &             - 0.5d0*(cmp(8,i)*cphi(9,i)+cmp(9,i)*cphi(8,i))
         vzz = vzz - cmp(4,i)*cphi(4,i) - 2.0d0*cmp(7,i)*cphi(7,i)
     &             - cmp(9,i)*cphi(9,i) - cmp(10,i)*cphi(10,i)
      end do
c
c     convert Cartesian induced dipoles to fractional coordinates
c

c      if (use_polar) then
c         do i = 1, 3
c            a(1,i) = dble(nfft1) * recip(i,1)
c            a(2,i) = dble(nfft2) * recip(i,2)
c            a(3,i) = dble(nfft3) * recip(i,3)
c         end do
c         do i = 1, npole
c            do j = 1, 3
c               fuind(j,i) = a(j,1)*uind(1,i) + a(j,2)*uind(2,i)
c     &                          + a(j,3)*uind(3,i)
c               fuinp(j,i) = a(j,1)*uinp(1,i) + a(j,2)*uinp(2,i)
c     &                          + a(j,3)*uinp(3,i)
c            end do
c         end do
c
c     assign PME grid and perform 3-D FFT forward transform
c
c         call grid_uind (fuind,fuinp)
c         call fftfront
c
c     account for the zeroth grid point for a finite system
c
c         if (.not. use_bounds) then
c            expterm = 0.5d0 * pi / xbox
c            struc2 = qgrid(1,1,1,1)**2 + qgrid(2,1,1,1)**2
c            e = 0.5d0 * expterm * struc2
c         end if
c
c     complete the transformation of the PME grid
c
c         do k = 1, nfft3
c            do j = 1, nfft2
c               do i = 1, nfft1
c                  term = qfac(i,j,k)
c                  qgrid(1,i,j,k) = term * qgrid(1,i,j,k)
c                  qgrid(2,i,j,k) = term * qgrid(2,i,j,k)
c               end do
c            end do
c         end do
c
c     perform 3-D FFT backward transform and get potential
c
c         call fftback
c         call fphi_uind (fphid,fphip,fphidp)
c         do i = 1, npole
c            do j = 1, 10
c               fphid(j,i) = electric * fphid(j,i)
c               fphip(j,i) = electric * fphip(j,i)
c            end do
c            do j = 1, 20
c               fphidp(j,i) = electric * fphidp(j,i)
c            end do
c         end do
c
c     increment the induced dipole energy and gradient
c
c         e = 0.0d0
c         do i = 1, npole
c            f1 = 0.0d0
c            f2 = 0.0d0
c            f3 = 0.0d0
c            do k = 1, 3
c               j1 = deriv1(k+1)
c               j2 = deriv2(k+1)
c               j3 = deriv3(k+1)
c               e = e + fuind(k,i)*fphi(k+1,i)
c               f1 = f1 + (fuind(k,i)+fuinp(k,i))*fphi(j1,i)
c     &                 + fuind(k,i)*fphip(j1,i)
c     &                 + fuinp(k,i)*fphid(j1,i)
c               f2 = f2 + (fuind(k,i)+fuinp(k,i))*fphi(j2,i)
c     &                 + fuind(k,i)*fphip(j2,i)
c     &                 + fuinp(k,i)*fphid(j2,i)
c               f3 = f3 + (fuind(k,i)+fuinp(k,i))*fphi(j3,i)
c     &                 + fuind(k,i)*fphip(j3,i)
c     &                 + fuinp(k,i)*fphid(j3,i)
c               if (poltyp .eq. 'DIRECT') then
c                  f1 = f1 - fuind(k,i)*fphip(j1,i)
c     &                    - fuinp(k,i)*fphid(j1,i)
c                  f2 = f2 - fuind(k,i)*fphip(j2,i)
c     &                    - fuinp(k,i)*fphid(j2,i)
c                  f3 = f3 - fuind(k,i)*fphip(j3,i)
c     &                    - fuinp(k,i)*fphid(j3,i)
c               end if
c            end do
c            do k = 1, 10
c               f1 = f1 + fmp(k,i)*fphidp(deriv1(k),i)
c               f2 = f2 + fmp(k,i)*fphidp(deriv2(k),i)
c               f3 = f3 + fmp(k,i)*fphidp(deriv3(k),i)
c            end do
c            f1 = 0.5d0 * dble(nfft1) * f1
c            f2 = 0.5d0 * dble(nfft2) * f2
c            f3 = 0.5d0 * dble(nfft3) * f3
c            frc(1,i) = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
c            frc(2,i) = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
c            frc(3,i) = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
c         end do
c         e = 0.5d0 * e
c         ep = ep + e
c         do i = 1, npole
c            ii = ipole(i)
c            dep(1,ii) = dep(1,ii) + frc(1,i)
c            dep(2,ii) = dep(2,ii) + frc(2,i)
c            dep(3,ii) = dep(3,ii) + frc(3,i)
c         end do
c
c     set the potential to be the induced dipole average
c
c         do i = 1, npole
c            do k = 1, 10
c               fphidp(k,i) = 0.5d0 * fphidp(k,i)
c            end do
c         end do
c         call fphi_to_cphi (fphidp,cphi)
c
c     distribute torques into the induced dipole gradient
c
c         do i = 1, npole
c            trq(1,i) = cmp(4,i)*cphi(3,i) - cmp(3,i)*cphi(4,i)
c     &                    + 2.0d0*(cmp(7,i)-cmp(6,i))*cphi(10,i)
c     &                    + cmp(9,i)*cphi(8,i) + cmp(10,i)*cphi(6,i)
c     &                    - cmp(8,i)*cphi(9,i) - cmp(10,i)*cphi(7,i)
c            trq(2,i) = cmp(2,i)*cphi(4,i) - cmp(4,i)*cphi(2,i)
c     &                    + 2.0d0*(cmp(5,i)-cmp(7,i))*cphi(9,i)
c     &                    + cmp(8,i)*cphi(10,i) + cmp(9,i)*cphi(7,i)
c     &                    - cmp(9,i)*cphi(5,i) - cmp(10,i)*cphi(8,i)
c            trq(3,i) = cmp(3,i)*cphi(2,i) - cmp(2,i)*cphi(3,i)
c     &                    + 2.0d0*(cmp(6,i)-cmp(5,i))*cphi(8,i)
c     &                    + cmp(8,i)*cphi(5,i) + cmp(10,i)*cphi(9,i)
c     &                    - cmp(8,i)*cphi(6,i) - cmp(9,i)*cphi(10,i)
c         end do
c         do i = 1, n
c            frc(1,i) = 0.0d0
c            frc(2,i) = 0.0d0
c            frc(3,i) = 0.0d0
c         end do
c         call torque2 (trq,frc)
c         do i = 1, n
c            dep(1,i) = dep(1,i) + frc(1,i)
c            dep(2,i) = dep(2,i) + frc(2,i)
c            dep(3,i) = dep(3,i) + frc(3,i)
c         end do
c
c     induced dipole contribution to the internal virial
c
c         do i = 1, npole
c            do j = 2, 4
c               cphim(j) = 0.0d0
c               cphid(j) = 0.0d0
c               cphip(j) = 0.0d0
c               do k = 2, 4
c                  cphim(j) = cphim(j) + ftc(j,k)*fphi(k,i)
c                  cphid(j) = cphid(j) + ftc(j,k)*fphid(k,i)
c                  cphip(j) = cphip(j) + ftc(j,k)*fphip(k,i)
c               end do
c            end do
c            vxx = vxx - cphi(2,i)*cmp(2,i)
c     &                - 0.5d0*(cphim(2)*(uind(1,i)+uinp(1,i))
c     &                        +cphid(2)*uinp(1,i)+cphip(2)*uind(1,i))
c            vyx = vyx - 0.5d0*(cphi(2,i)*cmp(3,i)+cphi(3,i)*cmp(2,i))
c     &                - 0.25d0*(cphim(2)*(uind(2,i)+uinp(2,i))
c     &                         +cphim(3)*(uind(1,i)+uinp(1,i))
c     &                         +cphid(2)*uinp(2,i)+cphip(2)*uind(2,i)
c     &                         +cphid(3)*uinp(1,i)+cphip(3)*uind(1,i))
c            vzx = vzx - 0.5d0*(cphi(2,i)*cmp(4,i)+cphi(4,i)*cmp(2,i))
c     &                - 0.25d0*(cphim(2)*(uind(3,i)+uinp(3,i))
c     &                         +cphim(4)*(uind(1,i)+uinp(1,i))
c     &                         +cphid(2)*uinp(3,i)+cphip(2)*uind(3,i)
c     &                         +cphid(4)*uinp(1,i)+cphip(4)*uind(1,i))
c            vyy = vyy - cphi(3,i)*cmp(3,i)
c     &                - 0.5d0*(cphim(3)*(uind(2,i)+uinp(2,i))
c     &                        +cphid(3)*uinp(2,i)+cphip(3)*uind(2,i))
c            vzy = vzy - 0.5d0*(cphi(3,i)*cmp(4,i)+cphi(4,i)*cmp(3,i))
c     &                - 0.25d0*(cphim(3)*(uind(3,i)+uinp(3,i))
c     &                         +cphim(4)*(uind(2,i)+uinp(2,i))
c     &                         +cphid(3)*uinp(3,i)+cphip(3)*uind(3,i)
c     &                         +cphid(4)*uinp(2,i)+cphip(4)*uind(2,i))
c            vzz = vzz - cphi(4,i)*cmp(4,i)
c     &                - 0.5d0*(cphim(4)*(uind(3,i)+uinp(3,i))
c     &                        +cphid(4)*uinp(3,i)+cphip(4)*uind(3,i))
c            vxx = vxx - 2.0d0*cmp(5,i)*cphi(5,i) - cmp(8,i)*cphi(8,i)
c     &                - cmp(9,i)*cphi(9,i)
c            vyx = vyx - (cmp(5,i)+cmp(6,i))*cphi(8,i)
c     &                - 0.5d0*(cmp(8,i)*(cphi(6,i)+cphi(5,i))
c     &                     +cmp(9,i)*cphi(10,i)+cmp(10,i)*cphi(9,i))
c            vzx = vzx - (cmp(5,i)+cmp(7,i))*cphi(9,i)
c     &                - 0.5d0*(cmp(9,i)*(cphi(5,i)+cphi(7,i))
c     &                     +cmp(8,i)*cphi(10,i)+cmp(10,i)*cphi(8,i))
c            vyy = vyy - 2.0d0*cmp(6,i)*cphi(6,i) - cmp(8,i)*cphi(8,i)
c     &                - cmp(10,i)*cphi(10,i)
c            vzy = vzy - (cmp(6,i)+cmp(7,i))*cphi(10,i)
c     &                - 0.5d0*(cmp(10,i)*(cphi(6,i)+cphi(7,i))
c     &                     +cmp(8,i)*cphi(9,i)+cmp(9,i)*cphi(8,i))
c            vzz = vzz - 2.0d0*cmp(7,i)*cphi(7,i) - cmp(9,i)*cphi(9,i)
c     &                - cmp(10,i)*cphi(10,i)
c            if (poltyp .eq. 'DIRECT') then
c               vxx = vxx + 0.5d0*(cphid(2)*uinp(1,i)+cphip(2)*uind(1,i))
c               vyx = vyx + 0.25d0*(cphid(2)*uinp(2,i)+cphip(2)*uind(2,i)
c     &                           +cphid(3)*uinp(1,i)+cphip(3)*uind(1,i))
c               vzx = vzx + 0.25d0*(cphid(2)*uinp(3,i)+cphip(2)*uind(3,i)
c     &                           +cphid(4)*uinp(1,i)+cphip(4)*uind(1,i))
c               vyy = vyy + 0.5d0*(cphid(3)*uinp(2,i)+cphip(3)*uind(2,i))
c               vzy = vzy + 0.25d0*(cphid(3)*uinp(3,i)+cphip(3)*uind(3,i)
c     &                           +cphid(4)*uinp(2,i)+cphip(4)*uind(2,i))
c               vzz = vzz + 0.5d0*(cphid(4)*uinp(3,i)+cphip(4)*uind(3,i))
c            end if
c         end do
c      end if

c
c     increment the internal virial tensor components
c
      vir(1,1) = vir(1,1) + vxx
      vir(2,1) = vir(2,1) + vyx
      vir(3,1) = vir(3,1) + vzx
      vir(1,2) = vir(1,2) + vyx
      vir(2,2) = vir(2,2) + vyy
      vir(3,2) = vir(3,2) + vzy
      vir(1,3) = vir(1,3) + vzx
      vir(2,3) = vir(2,3) + vzy
      vir(3,3) = vir(3,3) + vzz
c
c     perform deallocation of some local arrays
c
      deallocate (frc)
      deallocate (trq)
      deallocate (fuind)
      deallocate (fuinp)
      deallocate (cmp)
      deallocate (fmp)
      deallocate (fphi)
      deallocate (fphid)
      deallocate (fphip)
      deallocate (fphidp)
      deallocate (cphi)
      return
      end
