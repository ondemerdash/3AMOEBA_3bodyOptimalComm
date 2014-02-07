c
c     ################################################################
c     ##                                                            ##
c     ##               Subroutine ereal1c_3b_Perm                   ##
c     ##  Liam O-Suilleabhain, Omar Demerdash, Teresa Head-Gordon   ##
c     ##                 Spring 2013                                ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ereal1c_3b_Perm" evaluates the real space portion of the regular Ewald
c     summation energy and gradient due to atomic multipole interactions
c
c
      subroutine Innerloop_ereal1c_3b_Perm(atomind,emt,viremt,demt) 
      implicit none
      include 'sizes3b.i'
      include 'sizes.i'
      include 'atoms.i'
      include 'atoms3b.i'
      include 'bound.i'
      include 'boxes.i'
      include 'cell.i'
      include 'chgpot.i'
      include 'chgpot3b.i'
      include 'couple.i'
      include 'cutoff.i'
c      include 'deriv.i'
c      include 'energi.i'
c      include 'combo.i'
      include 'ewald.i'
      include 'math.i'
c      include 'molcul.i'
c      include 'molcul3b.i'
      include 'mplpot.i'
      include 'mpole3b.i'
      include 'mpole.i'
      include 'mpole3b_3.i'
      include 'mpole3b_2.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'shunt.i'
c      include 'virial.i'
      integer i,j,k,l1,l3
      integer ii,kk,jcell
      integer iax,iay,iaz
      integer kax,kay,kaz
      real*8 e,ei,f,bfac
      real*8 eintra,erfc
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 scale3,scale5
      real*8 scale7
      real*8 temp3,temp5,temp7
      real*8 dsc3,dsc5,dsc7
      real*8 psc3,psc5,psc7
      real*8 usc3,usc5
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 gfd,gfdr
      real*8 xr,yr,zr
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 xkx,ykx,zkx
      real*8 xky,yky,zky
      real*8 xkz,ykz,zkz
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 erl,erli
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 frcxi(3),frcxk(3)
      real*8 frcyi(3),frcyk(3)
      real*8 frczi(3),frczk(3)
      real*8 ci,di(3),qi(9)
      real*8 ck,dk(3),qk(9)
      real*8 fridmp(3),findmp(3)
      real*8 ftm2(3),ftm2i(3)
      real*8 ftm2r(3),ftm2ri(3)
      real*8 ttm2(3),ttm3(3)
      real*8 ttm2i(3),ttm3i(3)
      real*8 ttm2r(3),ttm3r(3)
      real*8 ttm2ri(3),ttm3ri(3)
      real*8 fdir(3),dixdk(3)
      real*8 dkxui(3),dixuk(3)
      real*8 dixukp(3),dkxuip(3)
      real*8 uixqkr(3),ukxqir(3)
      real*8 uixqkrp(3),ukxqirp(3)
      real*8 qiuk(3),qkui(3)
      real*8 qiukp(3),qkuip(3)
      real*8 rxqiuk(3),rxqkui(3)
      real*8 rxqiukp(3),rxqkuip(3)
      real*8 qidk(3),qkdi(3)
      real*8 qir(3),qkr(3)
      real*8 qiqkr(3),qkqir(3)
      real*8 qixqk(3),rxqir(3)
      real*8 dixr(3),dkxr(3)
      real*8 dixqkr(3),dkxqir(3)
      real*8 rxqkr(3),qkrxqir(3)
      real*8 rxqikr(3),rxqkir(3)
      real*8 rxqidk(3),rxqkdi(3)
      real*8 ddsc3(3),ddsc5(3)
      real*8 ddsc7(3)
      real*8 bn(0:5)
      real*8 sc(10),gl(0:8)
      real*8 sci(8),scip(8)
      real*8 gli(7),glip(7)
      real*8 gf(7),gfi(6)
      real*8 gfr(7),gfri(6)
      real*8 gti(6),gtri(6)
      real*8, allocatable :: mscale(:)
c      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
c      real*8, allocatable :: uscale(:)
      real*8 aewald
      character*6 mode
      integer atomind
      real*8 demt(3,*),viremt(3,3),emt
      external erfc
c
c
c     zero out the intramolecular portion of the Ewald energy
c
      aewald=aewaldPerm
      eintra = 0.0d0
      if (npole .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
c      allocate (pscale(n))
      allocate (dscale(n))
c      allocate (uscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         mscale(i) = 1.0d0
c         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
c         uscale(i) = 1.0d0
      end do

      emt=0.0d0
      do i=1,npole
         do j=1,3
           demt(j,i)=0.0d0
         end do 
      end do 
      do i=1,3
         do j=1,3
           viremt(i,j)=0.0d0
         end do
      end do 
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
c      mode = 'EWALD'
c      call switch (mode)
c
c     set the permanent multipole and induced dipole values
c

c      do i = 1, npole-1
         i=atomind
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         di(1) = rpole(2,i)
         di(2) = rpole(3,i)
         di(3) = rpole(4,i)
         qi(1) = rpole(5,i)
         qi(2) = rpole(6,i)
         qi(3) = rpole(7,i)
         qi(4) = rpole(8,i)
         qi(5) = rpole(9,i)
         qi(6) = rpole(10,i)
         qi(7) = rpole(11,i)
         qi(8) = rpole(12,i)
         qi(9) = rpole(13,i)
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
c            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
c            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
c            pscale(i14(j,ii)) = p4scale
c            do k = 1, np11(ii)
c                if (i14(j,ii) .eq. ip11(k,ii))
c     &            pscale(i14(j,ii)) = p4scale * p41scale
c            end do
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
c            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
c            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
c            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
c            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
c            uscale(ip14(j,ii)) = u4scale
         end do
         do k = i+1, npole
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dk(1) = rpole(2,k)
               dk(2) = rpole(3,k)
               dk(3) = rpole(4,k)
               qk(1) = rpole(5,k)
               qk(2) = rpole(6,k)
               qk(3) = rpole(7,k)
               qk(4) = rpole(8,k)
               qk(5) = rpole(9,k)
               qk(6) = rpole(10,k)
               qk(7) = rpole(11,k)
               qk(8) = rpole(12,k)
               qk(9) = rpole(13,k)
c
c     calculate the real space error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 5
                  bfac = dble(2*j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
c
c     apply Thole polarization damping to scale factors
c
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               rr11 = 9.0d0 * rr9 / r2
               scale3 = 1.0d0
               scale5 = 1.0d0
               scale7 = 1.0d0
               do j = 1, 3
                  ddsc3(j) = 0.0d0
                  ddsc5(j) = 0.0d0
                  ddsc7(j) = 0.0d0
               end do
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = 1.0d0 - expdamp
                     scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                     scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                       *expdamp
                     temp3 = -3.0d0 * damp * expdamp / r2
                     temp5 = -damp
                     temp7 = -0.2d0 - 0.6d0*damp
                     ddsc3(1) = temp3 * xr
                     ddsc3(2) = temp3 * yr
                     ddsc3(3) = temp3 * zr
                     ddsc5(1) = temp5 * ddsc3(1)
                     ddsc5(2) = temp5 * ddsc3(2)
                     ddsc5(3) = temp5 * ddsc3(3)
                     ddsc7(1) = temp7 * ddsc5(1)
                     ddsc7(2) = temp7 * ddsc5(2)
                     ddsc7(3) = temp7 * ddsc5(3)
                  end if
               end if
               dsc3 = 1.0d0 - scale3*dscale(kk)
               dsc5 = 1.0d0 - scale5*dscale(kk)
               dsc7 = 1.0d0 - scale7*dscale(kk)
c               psc3 = 1.0d0 - scale3*pscale(kk)
c               psc5 = 1.0d0 - scale5*pscale(kk)
c               psc7 = 1.0d0 - scale7*pscale(kk)
c               usc3 = 1.0d0 - scale3*uscale(kk)
c               usc5 = 1.0d0 - scale5*uscale(kk)

c
c     construct necessary auxiliary vectors
c
               dixdk(1) = di(2)*dk(3) - di(3)*dk(2)
               dixdk(2) = di(3)*dk(1) - di(1)*dk(3)
               dixdk(3) = di(1)*dk(2) - di(2)*dk(1)
               dixr(1) = di(2)*zr - di(3)*yr
               dixr(2) = di(3)*xr - di(1)*zr
               dixr(3) = di(1)*yr - di(2)*xr
               dkxr(1) = dk(2)*zr - dk(3)*yr
               dkxr(2) = dk(3)*xr - dk(1)*zr
               dkxr(3) = dk(1)*yr - dk(2)*xr
               qir(1) = qi(1)*xr + qi(4)*yr + qi(7)*zr
               qir(2) = qi(2)*xr + qi(5)*yr + qi(8)*zr
               qir(3) = qi(3)*xr + qi(6)*yr + qi(9)*zr
               qkr(1) = qk(1)*xr + qk(4)*yr + qk(7)*zr
               qkr(2) = qk(2)*xr + qk(5)*yr + qk(8)*zr
               qkr(3) = qk(3)*xr + qk(6)*yr + qk(9)*zr
               qiqkr(1) = qi(1)*qkr(1) + qi(4)*qkr(2) + qi(7)*qkr(3)
               qiqkr(2) = qi(2)*qkr(1) + qi(5)*qkr(2) + qi(8)*qkr(3)
               qiqkr(3) = qi(3)*qkr(1) + qi(6)*qkr(2) + qi(9)*qkr(3)
               qkqir(1) = qk(1)*qir(1) + qk(4)*qir(2) + qk(7)*qir(3)
               qkqir(2) = qk(2)*qir(1) + qk(5)*qir(2) + qk(8)*qir(3)
               qkqir(3) = qk(3)*qir(1) + qk(6)*qir(2) + qk(9)*qir(3)
               qixqk(1) = qi(2)*qk(3) + qi(5)*qk(6) + qi(8)*qk(9)
     &                       - qi(3)*qk(2) - qi(6)*qk(5) - qi(9)*qk(8)
               qixqk(2) = qi(3)*qk(1) + qi(6)*qk(4) + qi(9)*qk(7)
     &                       - qi(1)*qk(3) - qi(4)*qk(6) - qi(7)*qk(9)
               qixqk(3) = qi(1)*qk(2) + qi(4)*qk(5) + qi(7)*qk(8)
     &                       - qi(2)*qk(1) - qi(5)*qk(4) - qi(8)*qk(7)
               rxqir(1) = yr*qir(3) - zr*qir(2)
               rxqir(2) = zr*qir(1) - xr*qir(3)
               rxqir(3) = xr*qir(2) - yr*qir(1)
               rxqkr(1) = yr*qkr(3) - zr*qkr(2)
               rxqkr(2) = zr*qkr(1) - xr*qkr(3)
               rxqkr(3) = xr*qkr(2) - yr*qkr(1)
               rxqikr(1) = yr*qiqkr(3) - zr*qiqkr(2)
               rxqikr(2) = zr*qiqkr(1) - xr*qiqkr(3)
               rxqikr(3) = xr*qiqkr(2) - yr*qiqkr(1)
               rxqkir(1) = yr*qkqir(3) - zr*qkqir(2)
               rxqkir(2) = zr*qkqir(1) - xr*qkqir(3)
               rxqkir(3) = xr*qkqir(2) - yr*qkqir(1)
               qkrxqir(1) = qkr(2)*qir(3) - qkr(3)*qir(2)
               qkrxqir(2) = qkr(3)*qir(1) - qkr(1)*qir(3)
               qkrxqir(3) = qkr(1)*qir(2) - qkr(2)*qir(1)
               qidk(1) = qi(1)*dk(1) + qi(4)*dk(2) + qi(7)*dk(3)
               qidk(2) = qi(2)*dk(1) + qi(5)*dk(2) + qi(8)*dk(3)
               qidk(3) = qi(3)*dk(1) + qi(6)*dk(2) + qi(9)*dk(3)
               qkdi(1) = qk(1)*di(1) + qk(4)*di(2) + qk(7)*di(3)
               qkdi(2) = qk(2)*di(1) + qk(5)*di(2) + qk(8)*di(3)
               qkdi(3) = qk(3)*di(1) + qk(6)*di(2) + qk(9)*di(3)
               dixqkr(1) = di(2)*qkr(3) - di(3)*qkr(2)
               dixqkr(2) = di(3)*qkr(1) - di(1)*qkr(3)
               dixqkr(3) = di(1)*qkr(2) - di(2)*qkr(1)
               dkxqir(1) = dk(2)*qir(3) - dk(3)*qir(2)
               dkxqir(2) = dk(3)*qir(1) - dk(1)*qir(3)
               dkxqir(3) = dk(1)*qir(2) - dk(2)*qir(1)
               rxqidk(1) = yr*qidk(3) - zr*qidk(2)
               rxqidk(2) = zr*qidk(1) - xr*qidk(3)
               rxqidk(3) = xr*qidk(2) - yr*qidk(1)
               rxqkdi(1) = yr*qkdi(3) - zr*qkdi(2)
               rxqkdi(2) = zr*qkdi(1) - xr*qkdi(3)
               rxqkdi(3) = xr*qkdi(2) - yr*qkdi(1)
c
c     calculate the scalar products for permanent components
c
               sc(2) = di(1)*dk(1) + di(2)*dk(2) + di(3)*dk(3)
               sc(3) = di(1)*xr + di(2)*yr + di(3)*zr
               sc(4) = dk(1)*xr + dk(2)*yr + dk(3)*zr
               sc(5) = qir(1)*xr + qir(2)*yr + qir(3)*zr
               sc(6) = qkr(1)*xr + qkr(2)*yr + qkr(3)*zr
               sc(7) = qir(1)*dk(1) + qir(2)*dk(2) + qir(3)*dk(3)
               sc(8) = qkr(1)*di(1) + qkr(2)*di(2) + qkr(3)*di(3)
               sc(9) = qir(1)*qkr(1) + qir(2)*qkr(2) + qir(3)*qkr(3)
               sc(10) = qi(1)*qk(1) + qi(2)*qk(2) + qi(3)*qk(3)
     &                     + qi(4)*qk(4) + qi(5)*qk(5) + qi(6)*qk(6)
     &                     + qi(7)*qk(7) + qi(8)*qk(8) + qi(9)*qk(9)
c
c     calculate the gl functions for permanent components
c
               gl(0) = ci*ck
               gl(1) = ck*sc(3) - ci*sc(4)
               gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
               gl(3) = sc(3)*sc(6) - sc(4)*sc(5)
               gl(4) = sc(5)*sc(6)
               gl(5) = -4.0d0 * sc(9)
               gl(6) = sc(2)
               gl(7) = 2.0d0 * (sc(7)-sc(8))
               gl(8) = 2.0d0 * sc(10)
c
c     compute the energy contributions for this interaction
c
               e = bn(0)*gl(0) + bn(1)*(gl(1)+gl(6))
     &                + bn(2)*(gl(2)+gl(7)+gl(8))
     &                + bn(3)*(gl(3)+gl(5)) + bn(4)*gl(4)
c
c     get the real energy without any screening function
c
               erl = rr1*gl(0) + rr3*(gl(1)+gl(6))
     &                  + rr5*(gl(2)+gl(7)+gl(8))
     &                  + rr7*(gl(3)+gl(5)) + rr9*gl(4)
               e = e - (1.0d0-mscale(kk))*erl
               e = f * e
c               em = em + e
               emt = emt + e 
c
c     increment the total intramolecular energy; assumes
c     intramolecular distances are less than half of cell
c     length and less than the ewald cutoff
c
               if (molcule(ii) .eq. molcule(kk)) then
                 eintra = eintra + mscale(kk)*erl*f
               end if
c
c     intermediate variables for permanent force terms
c
               gf(1) = bn(1)*gl(0) + bn(2)*(gl(1)+gl(6))
     &                    + bn(3)*(gl(2)+gl(7)+gl(8))
     &                    + bn(4)*(gl(3)+gl(5)) + bn(5)*gl(4)
               gf(2) = -ck*bn(1) + sc(4)*bn(2) - sc(6)*bn(3)
               gf(3) = ci*bn(1) + sc(3)*bn(2) + sc(5)*bn(3)
               gf(4) = 2.0d0 * bn(2)
               gf(5) = 2.0d0 * (-ck*bn(2)+sc(4)*bn(3)-sc(6)*bn(4))
               gf(6) = 2.0d0 * (-ci*bn(2)-sc(3)*bn(3)-sc(5)*bn(4))
               gf(7) = 4.0d0 * bn(3)
               gfr(1) = rr3*gl(0) + rr5*(gl(1)+gl(6))
     &                     + rr7*(gl(2)+gl(7)+gl(8))
     &                     + rr9*(gl(3)+gl(5)) + rr11*gl(4)
               gfr(2) = -ck*rr3 + sc(4)*rr5 - sc(6)*rr7
               gfr(3) = ci*rr3 + sc(3)*rr5 + sc(5)*rr7
               gfr(4) = 2.0d0 * rr5
               gfr(5) = 2.0d0 * (-ck*rr5+sc(4)*rr7-sc(6)*rr9)
               gfr(6) = 2.0d0 * (-ci*rr5-sc(3)*rr7-sc(5)*rr9)
               gfr(7) = 4.0d0 * rr7
c
c     get the permanent force with screening
c
               ftm2(1) = gf(1)*xr + gf(2)*di(1) + gf(3)*dk(1)
     &                      + gf(4)*(qkdi(1)-qidk(1)) + gf(5)*qir(1)
     &                      + gf(6)*qkr(1) + gf(7)*(qiqkr(1)+qkqir(1))
               ftm2(2) = gf(1)*yr + gf(2)*di(2) + gf(3)*dk(2)
     &                      + gf(4)*(qkdi(2)-qidk(2)) + gf(5)*qir(2)
     &                      + gf(6)*qkr(2) + gf(7)*(qiqkr(2)+qkqir(2))
               ftm2(3) = gf(1)*zr + gf(2)*di(3) + gf(3)*dk(3)
     &                      + gf(4)*(qkdi(3)-qidk(3)) + gf(5)*qir(3)
     &                      + gf(6)*qkr(3) + gf(7)*(qiqkr(3)+qkqir(3))
c
c     get the permanent force without screening
c
               ftm2r(1) = gfr(1)*xr + gfr(2)*di(1) + gfr(3)*dk(1)
     &                      + gfr(4)*(qkdi(1)-qidk(1)) + gfr(5)*qir(1)
     &                      + gfr(6)*qkr(1) + gfr(7)*(qiqkr(1)+qkqir(1))
               ftm2r(2) = gfr(1)*yr + gfr(2)*di(2) + gfr(3)*dk(2)
     &                      + gfr(4)*(qkdi(2)-qidk(2)) + gfr(5)*qir(2)
     &                      + gfr(6)*qkr(2) + gfr(7)*(qiqkr(2)+qkqir(2))
               ftm2r(3) = gfr(1)*zr + gfr(2)*di(3) + gfr(3)*dk(3)
     &                      + gfr(4)*(qkdi(3)-qidk(3)) + gfr(5)*qir(3)
     &                      + gfr(6)*qkr(3) + gfr(7)*(qiqkr(3)+qkqir(3))

c
c     account for partially excluded induced interactions
c
c               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(kk)
c     &                                  +(glip(1)+glip(6))*dscale(kk))
c               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(kk)
c     &                                  +(glip(2)+glip(7))*dscale(kk))
c               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(kk)
c     &                                  +glip(3)*dscale(kk))
               fridmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
     &                        + temp7*ddsc7(1)
               fridmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
     &                        + temp7*ddsc7(2)
               fridmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
     &                        + temp7*ddsc7(3)
c
c     get the permanent torque with screening
c
               ttm2(1) = -bn(1)*dixdk(1) + gf(2)*dixr(1)
     &           + gf(4)*(dixqkr(1)+dkxqir(1)+rxqidk(1)-2.0d0*qixqk(1))
     &           - gf(5)*rxqir(1) - gf(7)*(rxqikr(1)+qkrxqir(1))
               ttm2(2) = -bn(1)*dixdk(2) + gf(2)*dixr(2)
     &           + gf(4)*(dixqkr(2)+dkxqir(2)+rxqidk(2)-2.0d0*qixqk(2))
     &           - gf(5)*rxqir(2) - gf(7)*(rxqikr(2)+qkrxqir(2))
               ttm2(3) = -bn(1)*dixdk(3) + gf(2)*dixr(3)
     &           + gf(4)*(dixqkr(3)+dkxqir(3)+rxqidk(3)-2.0d0*qixqk(3))
     &           - gf(5)*rxqir(3) - gf(7)*(rxqikr(3)+qkrxqir(3))
               ttm3(1) = bn(1)*dixdk(1) + gf(3)*dkxr(1)
     &           - gf(4)*(dixqkr(1)+dkxqir(1)+rxqkdi(1)-2.0d0*qixqk(1))
     &           - gf(6)*rxqkr(1) - gf(7)*(rxqkir(1)-qkrxqir(1))
               ttm3(2) = bn(1)*dixdk(2) + gf(3)*dkxr(2)
     &           - gf(4)*(dixqkr(2)+dkxqir(2)+rxqkdi(2)-2.0d0*qixqk(2))
     &           - gf(6)*rxqkr(2) - gf(7)*(rxqkir(2)-qkrxqir(2))
               ttm3(3) = bn(1)*dixdk(3) + gf(3)*dkxr(3)
     &           - gf(4)*(dixqkr(3)+dkxqir(3)+rxqkdi(3)-2.0d0*qixqk(3))
     &           - gf(6)*rxqkr(3) - gf(7)*(rxqkir(3)-qkrxqir(3))
c
c     get the permanent torque without screening
c
               ttm2r(1) = -rr3*dixdk(1) + gfr(2)*dixr(1)-gfr(5)*rxqir(1)
     &           + gfr(4)*(dixqkr(1)+dkxqir(1)+rxqidk(1)-2.0d0*qixqk(1))
     &           - gfr(7)*(rxqikr(1)+qkrxqir(1))
               ttm2r(2) = -rr3*dixdk(2) + gfr(2)*dixr(2)-gfr(5)*rxqir(2)
     &           + gfr(4)*(dixqkr(2)+dkxqir(2)+rxqidk(2)-2.0d0*qixqk(2))
     &           - gfr(7)*(rxqikr(2)+qkrxqir(2))
               ttm2r(3) = -rr3*dixdk(3) + gfr(2)*dixr(3)-gfr(5)*rxqir(3)
     &           + gfr(4)*(dixqkr(3)+dkxqir(3)+rxqidk(3)-2.0d0*qixqk(3))
     &           - gfr(7)*(rxqikr(3)+qkrxqir(3))
               ttm3r(1) = rr3*dixdk(1) + gfr(3)*dkxr(1) -gfr(6)*rxqkr(1)
     &           - gfr(4)*(dixqkr(1)+dkxqir(1)+rxqkdi(1)-2.0d0*qixqk(1))
     &           - gfr(7)*(rxqkir(1)-qkrxqir(1))
               ttm3r(2) = rr3*dixdk(2) + gfr(3)*dkxr(2) -gfr(6)*rxqkr(2)
     &           - gfr(4)*(dixqkr(2)+dkxqir(2)+rxqkdi(2)-2.0d0*qixqk(2))
     &           - gfr(7)*(rxqkir(2)-qkrxqir(2))
               ttm3r(3) = rr3*dixdk(3) + gfr(3)*dkxr(3) -gfr(6)*rxqkr(3)
     &           - gfr(4)*(dixqkr(3)+dkxqir(3)+rxqkdi(3)-2.0d0*qixqk(3))
     &           - gfr(7)*(rxqkir(3)-qkrxqir(3))
c
c     handle the case where scaling is used
c
               do j = 1, 3
                  ftm2(j) = f * (ftm2(j)-(1.0d0-mscale(kk))*ftm2r(j))
                  ttm2(j) = f * (ttm2(j)-(1.0d0-mscale(kk))*ttm2r(j))
                  ttm3(j) = f * (ttm3(j)-(1.0d0-mscale(kk))*ttm3r(j))
               end do
c
c     increment gradient due to force and torque on first site
c

c               dem(1,ii) = dem(1,ii) + ftm2(1)
c               dem(2,ii) = dem(2,ii) + ftm2(2)
c               dem(3,ii) = dem(3,ii) + ftm2(3)

               demt(1,ii) = demt(1,ii) + ftm2(1)
               demt(2,ii) = demt(2,ii) + ftm2(2)
               demt(3,ii) = demt(3,ii) + ftm2(3)
               call torque_3b_Perm_new (demt,
     &         i,ttm2,ttm2i,frcxi,frcyi,frczi)
c
c     increment gradient due to force and torque on second site
c
c               dem(1,kk) = dem(1,kk) - ftm2(1)
c               dem(2,kk) = dem(2,kk) - ftm2(2)
c               dem(3,kk) = dem(3,kk) - ftm2(3)

               demt(1,kk) = demt(1,kk) - ftm2(1)
               demt(2,kk) = demt(2,kk) - ftm2(2)
               demt(3,kk) = demt(3,kk) - ftm2(3)
               call torque_3b_Perm_new (demt,
     &          k,ttm3,ttm3i,frcxk,frcyk,frczk)

c
c     increment the internal virial tensor components
c
               iaz = zaxis(i)
               iax = xaxis(i)
               iay = yaxis(i)
               kaz = zaxis(k)
               kax = xaxis(k)
               kay = yaxis(k)
               if (iaz .eq. 0)  iaz = ii
               if (iax .eq. 0)  iax = ii
               if (iay .eq. 0)  iay = ii
               if (kaz .eq. 0)  kaz = kk
               if (kax .eq. 0)  kax = kk
               if (kay .eq. 0)  kay = kk
               xiz = x(iaz) - x(ii)
               yiz = y(iaz) - y(ii)
               ziz = z(iaz) - z(ii)
               xix = x(iax) - x(ii)
               yix = y(iax) - y(ii)
               zix = z(iax) - z(ii)
               xiy = x(iay) - x(ii)
               yiy = y(iay) - y(ii)
               ziy = z(iay) - z(ii)
               xkz = x(kaz) - x(kk)
               ykz = y(kaz) - y(kk)
               zkz = z(kaz) - z(kk)
               xkx = x(kax) - x(kk)
               ykx = y(kax) - y(kk)
               zkx = z(kax) - z(kk)
               xky = x(kay) - x(kk)
               yky = y(kay) - y(kk)
               zky = z(kay) - z(kk)
               vxx = -xr*(ftm2(1)) + xix*frcxi(1)
     &                  + xiy*frcyi(1) + xiz*frczi(1) + xkx*frcxk(1)
     &                  + xky*frcyk(1) + xkz*frczk(1)
               vyx = -yr*(ftm2(1)) + yix*frcxi(1)
     &                  + yiy*frcyi(1) + yiz*frczi(1) + ykx*frcxk(1)
     &                  + yky*frcyk(1) + ykz*frczk(1)
               vzx = -zr*(ftm2(1)) + zix*frcxi(1)
     &                  + ziy*frcyi(1) + ziz*frczi(1) + zkx*frcxk(1)
     &                  + zky*frcyk(1) + zkz*frczk(1)
               vyy = -yr*(ftm2(2)) + yix*frcxi(2)
     &                  + yiy*frcyi(2) + yiz*frczi(2) + ykx*frcxk(2)
     &                  + yky*frcyk(2) + ykz*frczk(2)
               vzy = -zr*(ftm2(2)) + zix*frcxi(2)
     &                  + ziy*frcyi(2) + ziz*frczi(2) + zkx*frcxk(2)
     &                  + zky*frcyk(2) + zkz*frczk(2)
               vzz = -zr*(ftm2(3)) + zix*frcxi(3)
     &                  + ziy*frcyi(3) + ziz*frczi(3) + zkx*frcxk(3)
     &                  + zky*frcyk(3) + zkz*frczk(3)
c               vir(1,1) = vir(1,1) + vxx
c               vir(2,1) = vir(2,1) + vyx
c               vir(3,1) = vir(3,1) + vzx
c               vir(1,2) = vir(1,2) + vyx
c               vir(2,2) = vir(2,2) + vyy
c               vir(3,2) = vir(3,2) + vzy
c               vir(1,3) = vir(1,3) + vzx
c               vir(2,3) = vir(2,3) + vzy
c               vir(3,3) = vir(3,3) + vzz

               viremt(1,1) = viremt(1,1) + vxx
               viremt(2,1) = viremt(2,1) + vyx
               viremt(3,1) = viremt(3,1) + vzx
               viremt(1,2) = viremt(1,2) + vyx
               viremt(2,2) = viremt(2,2) + vyy
               viremt(3,2) = viremt(3,2) + vzy
               viremt(1,3) = viremt(1,3) + vzx
               viremt(2,3) = viremt(2,3) + vzy
               viremt(3,3) = viremt(3,3) + vzz
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
c            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
c            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
c            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
c            pscale(i15(j,ii)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
c            uscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
c            uscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
c            uscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
c            uscale(ip14(j,ii)) = 1.0d0
         end do
c      end do
c
c     set the permanent multipole and induced dipole values
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
c      deallocate (pscale)
      deallocate (dscale)
c      deallocate (uscale)
      return
      end
