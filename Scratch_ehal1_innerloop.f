c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ehal1a  --  double loop buffer 14-7 vdw derivs  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ehal1a" calculates the buffered 14-7 van der Waals energy and
c     its first derivatives with respect to Cartesian coordinates
c     using a pairwise double loop
c
c
c      subroutine ehal1a
      subroutine Innerloop_ehal1a(atomind,evt,virevt,devt)
      implicit none
      include 'sizes3b.i'
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'atoms3b.i'
      include 'bound.i'
      include 'cell.i'
      include 'couple.i'
c      include 'deriv.i'
c      include 'energi.i'
      include 'group.i'
      include 'group3b.i'
      include 'inter.i'
      include 'molcul.i'
      include 'molcul3b.i'
      include 'shunt.i'
      include 'usage.i'
      include 'vdw.i'
      include 'vdwpot.i'
c      include 'virial.i'
      integer i,j,k
      integer ii,iv,it
      integer kk,kv,kt
      integer, allocatable :: iv14(:)
      real*8 e,de,eps,rdn
      real*8 fgrp,rv,rv7
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 redi,rediv
      real*8 redk,redkv
      real*8 dedx,dedy,dedz
      real*8 rho,tau,tau7
      real*8 dtau,gtau
      real*8 taper,dtaper
      real*8 rik,rik2,rik3
      real*8 rik4,rik5
      real*8 rik6,rik7
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
      real*8, allocatable :: vscale(:)
      logical proceed,usei
      character*6 mode
      integer atomind
      real*8 evt,virevt(3,3),devt(3,*)
c
c
c     zero out the van der Waals energy and first derivatives
c
      evt = 0.0d0
      do i = 1, n
         devt(1,i) = 0.0d0
         devt(2,i) = 0.0d0
         devt(3,i) = 0.0d0
      end do
      do i=1,3
         do j=1,3
           virevt(i,j)=0.0d0
         end do 
      end do 
c
c     perform dynamic allocation of some local arrays
c
      allocate (iv14(n))
      allocate (xred(n))
      allocate (yred(n))
      allocate (zred(n))
      allocate (vscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         vscale(i) = 1.0d0
         iv14(i) = 0
      end do
c
c     set the coefficients for the switching function
c

c      mode = 'VDW'
c      call switch (mode)

c
c     apply any reduction factor to the atomic coordinates
c
      do k = 1, nvdw
         i = ivdw(k)
         iv = ired(i)
         rdn = kred(i)
         xred(i) = rdn*(x(i)-x(iv)) + x(iv)
         yred(i) = rdn*(y(i)-y(iv)) + y(iv)
         zred(i) = rdn*(z(i)-z(iv)) + z(iv)
      end do
c
c     find van der Waals energy and derivatives via double loop
c

c      do ii = 1, nvdw-1
         ii=atomind
         i = ivdw(ii)
         iv = ired(i)
         redi = kred(i)
         rediv = 1.0d0 - redi
         it = jvdw(i)
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         usei = (use(i) .or. use(iv))
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(i)
            vscale(i12(j,i)) = v2scale
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = v3scale
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = v4scale
            iv14(i14(j,i)) = i
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = v5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii+1, nvdw
            k = ivdw(kk)
            kv = ired(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kv))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               kt = jvdw(k)
               xr = xi - xred(k)
               yr = yi - yred(k)
               zr = zi - zred(k)
               call image (xr,yr,zr)
               rik2 = xr*xr + yr*yr + zr*zr
c
c     check for an interaction distance less than the cutoff
c
               if (rik2 .le. off2) then
                  rv = radmin(kt,it)
                  eps = epsilon(kt,it)
                  if (iv14(k) .eq. i) then
                     rv = radmin4(kt,it)
                     eps = epsilon4(kt,it)
                  end if
                  eps = eps * vscale(k)
                  rik = sqrt(rik2)
                  rv7 = rv**7
                  rik6 = rik2**3
                  rik7 = rik6 * rik
                  rho = rik7 + ghal*rv7
                  tau = (dhal+1.0d0) / (rik + dhal*rv)
                  tau7 = tau**7
                  dtau = tau / (dhal+1.0d0)
                  gtau = eps*tau7*rik6*(ghal+1.0d0)*(rv7/rho)**2
                  e = eps*tau7*rv7*((ghal+1.0d0)*rv7/rho-2.0d0)
                  de = -7.0d0 * (dtau*e+gtau)
c
c     use energy switching if near the cutoff distance
c
                  if (rik2 .gt. cut2) then
                     rik3 = rik2 * rik
                     rik4 = rik2 * rik2
                     rik5 = rik2 * rik3
                     taper = c5*rik5 + c4*rik4 + c3*rik3
     &                          + c2*rik2 + c1*rik + c0
                     dtaper = 5.0d0*c5*rik4 + 4.0d0*c4*rik3
     &                           + 3.0d0*c3*rik2 + 2.0d0*c2*rik + c1
                     de = e*dtaper + de*taper
                     e = e * taper
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                  end if
c
c     find the chain rule terms for derivative components
c
                  de = de / rik
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
c
c     increment the total van der Waals energy and derivatives
c
c                  ev = ev + e
                  evt = evt + e

                  if (i .eq. iv) then
c                     dev(1,i) = dev(1,i) + dedx
c                     dev(2,i) = dev(2,i) + dedy
c                     dev(3,i) = dev(3,i) + dedz
                     devt(1,i) = devt(1,i) + dedx
                     devt(2,i) = devt(2,i) + dedy
                     devt(3,i) = devt(3,i) + dedz

                  else
c                     dev(1,i) = dev(1,i) + dedx*redi
c                     dev(2,i) = dev(2,i) + dedy*redi
c                     dev(3,i) = dev(3,i) + dedz*redi
c                     dev(1,iv) = dev(1,iv) + dedx*rediv
c                     dev(2,iv) = dev(2,iv) + dedy*rediv
c                     dev(3,iv) = dev(3,iv) + dedz*rediv

                     devt(1,i) = devt(1,i) + dedx*redi
                     devt(2,i) = devt(2,i) + dedy*redi
                     devt(3,i) = devt(3,i) + dedz*redi
                     devt(1,iv) = devt(1,iv) + dedx*rediv
                     devt(2,iv) = devt(2,iv) + dedy*rediv
                     devt(3,iv) = devt(3,iv) + dedz*rediv

                  end if
                  if (k .eq. kv) then
c                     dev(1,k) = dev(1,k) - dedx
c                     dev(2,k) = dev(2,k) - dedy
c                     dev(3,k) = dev(3,k) - dedz

                     devt(1,k) = devt(1,k) - dedx
                     devt(2,k) = devt(2,k) - dedy
                     devt(3,k) = devt(3,k) - dedz

                  else
                     redk = kred(k)
                     redkv = 1.0d0 - redk
c                     dev(1,k) = dev(1,k) - dedx*redk
c                     dev(2,k) = dev(2,k) - dedy*redk
c                     dev(3,k) = dev(3,k) - dedz*redk
c                     dev(1,kv) = dev(1,kv) - dedx*redkv
c                     dev(2,kv) = dev(2,kv) - dedy*redkv
c                     dev(3,kv) = dev(3,kv) - dedz*redkv

                     devt(1,k) = devt(1,k) - dedx*redk
                     devt(2,k) = devt(2,k) - dedy*redk
                     devt(3,k) = devt(3,k) - dedz*redk
                     devt(1,kv) = devt(1,kv) - dedx*redkv
                     devt(2,kv) = devt(2,kv) - dedy*redkv
                     devt(3,kv) = devt(3,kv) - dedz*redkv
                  end if
c
c     increment the internal virial tensor components
c
                  vxx = xr * dedx
                  vyx = yr * dedx
                  vzx = zr * dedx
                  vyy = yr * dedy
                  vzy = zr * dedy
                  vzz = zr * dedz
c                  vir(1,1) = vir(1,1) + vxx
c                  vir(2,1) = vir(2,1) + vyx
c                  vir(3,1) = vir(3,1) + vzx
c                  vir(1,2) = vir(1,2) + vyx
c                  vir(2,2) = vir(2,2) + vyy
c                  vir(3,2) = vir(3,2) + vzy
c                  vir(1,3) = vir(1,3) + vzx
c                  vir(2,3) = vir(2,3) + vzy
c                  vir(3,3) = vir(3,3) + vzz

                  virevt(1,1) = virevt(1,1) + vxx
                  virevt(2,1) = virevt(2,1) + vyx
                  virevt(3,1) = virevt(3,1) + vzx
                  virevt(1,2) = virevt(1,2) + vyx
                  virevt(2,2) = virevt(2,2) + vyy
                  virevt(3,2) = virevt(3,2) + vzy
                  virevt(1,3) = virevt(1,3) + vzx
                  virevt(2,3) = virevt(2,3) + vzy
                  virevt(3,3) = virevt(3,3) + vzz

c
c     increment the total intermolecular energy
c
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(i)
            vscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            vscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            vscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            vscale(i15(j,i)) = 1.0d0
         end do
c      end do

      deallocate (iv14)
      deallocate (xred)
      deallocate (yred)
      deallocate (zred)
      deallocate (vscale)
      return
      end

