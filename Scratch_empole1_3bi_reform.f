c
c     ################################################################
c     ##                                                            ##
c     ##              Subroutine empole1c_3b_Perm                   ##
c     ##  Liam O-Suilleabhain, Omar Demerdash, Teresa Head-Gordon   ##
c     ##                 Spring 2013                                ##
c     ##                                                            ##
c     ################################################################
c
c
c     "empole1c_3b_Perm" calculates the multipole energy and derivatives 
c     with respect to Cartesian coordinates using regular Ewald 
c
c
c      subroutine empole1c_3b_Perm(field,fieldp)
      subroutine empole1c_3b_Perm_PME
      implicit none
      include 'sizes3b.i'
      include 'sizes.i'
      include 'atoms.i'
      include 'atoms3b.i'
      include 'boxes.i'
      include 'combo.i'
      include 'chgpot.i'
      include 'chgpot3b.i'
      include 'deriv.i'
      include 'energi.i'
      include 'ewald.i'
      include 'inter.i'
      include 'math.i'
      include 'mpole3b.i'
      include 'mpole.i'
      include 'mpole3b_3.i'
      include 'mpole3b_2.i'
      include 'polar.i'
      include 'potent.i'
      include 'polpot.i'
      include 'virial.i'
      integer i,j,ii
      real*8 e,ei,eintra
      real*8 f,term,fterm
      real*8 cii,dii,qii,uii
      real*8 xd,yd,zd
      real*8 xu,yu,zu
      real*8 xup,yup,zup
      real*8 xq,yq,zq
      real*8 xv,yv,zv,vterm
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 xdfield,xufield
      real*8 ydfield,yufield
      real*8 zdfield,zufield
      real*8 trq(3),trqi(3)
      real*8 frcx(3),frcy(3),frcz(3)
      real*8 aewald

      aewald=aewaldPerm
c      real*8 field(3,npole)
c      real*8 fieldp(3,npole)
c
c     zero out multipole energy and derivatives
c
      em = 0.0d0
      do i = 1, n
        do j = 1, 3
          dem(j,i) = 0.0d0
        end do
      end do
      if (.not.use_mpole)return
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     get the field at each atom
c

c      call field_3b_Perm(field,fieldp)

c
c     compute the reciprocal space part of the Ewald summation
c

c      call erecip1_3b_Perm
      call emrecip1_3b_Perm

      print*,"After erecip1_3b_Perm"
c
c     compute the real space part of the Ewald summation
c
      call ereal1c_3b_Perm(eintra)
      print*,"After ereal1c_3b_Perm"
c
c     compute the Ewald self-energy term over all the atoms
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      do i = 1, npole
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         cii = ci*ci
         dii = dix*dix + diy*diy + diz*diz
         qii = qixx*qixx + qiyy*qiyy + qizz*qizz
     &            + 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
         e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
         em = em + e
      end do
c
c     compute the cell dipole boundary correction term
c

      if (boundary .eq. 'VACUUM') then
        trqi(1) = 0.0d0
        trqi(2) = 0.0d0
        trqi(3) = 0.0d0
        xd = 0.0d0
        yd = 0.0d0
        zd = 0.0d0
        xu = 0.0d0
        yu = 0.0d0
        zu = 0.0d0
        xup = 0.0d0
        yup = 0.0d0
        zup = 0.0d0
        do i = 1, npole
          ii = ipole(i)
          xd = xd + rpole(2,i) + rpole(1,i)*x(ii)
          yd = yd + rpole(3,i) + rpole(1,i)*y(ii)
          zd = zd + rpole(4,i) + rpole(1,i)*z(ii)
        end do
        term = (2.0d0/3.0d0) * f * (pi/volbox)
        em = em + term*(xd*xd+yd*yd+zd*zd)
        do i = 1, npole
          ii = ipole(i)
          dem(1,ii) = dem(1,ii) + 2.0d0*term*rpole(1,i)*xd
          dem(2,ii) = dem(2,ii) + 2.0d0*term*rpole(1,i)*yd
          dem(3,ii) = dem(3,ii) + 2.0d0*term*rpole(1,i)*zd
        end do
        xdfield = -2.0d0 * term * xd
        ydfield = -2.0d0 * term * yd
        zdfield = -2.0d0 * term * zd
        do i = 1, npole
          trq(1) = rpole(3,i)*zdfield - rpole(4,i)*ydfield
          trq(2) = rpole(4,i)*xdfield - rpole(2,i)*zdfield
          trq(3) = rpole(2,i)*ydfield - rpole(3,i)*xdfield
          call torque_3b_Perm (i,trq,trqi,frcx,frcy,frcz)
        end do
c
c     boundary correction to virial due to overall cell dipole
c
        xd = 0.0d0
        yd = 0.0d0
        zd = 0.0d0
        xq = 0.0d0
        yq = 0.0d0
        zq = 0.0d0
        do i = 1, npole
          ii = ipole(i)
          xd = xd + rpole(2,i)
          yd = yd + rpole(3,i)
          zd = zd + rpole(4,i)
          xq = xq + rpole(1,i)*x(ii)
          yq = yq + rpole(1,i)*y(ii)
          zq = zq + rpole(1,i)*z(ii)
        end do
        xv = xq * (xd+0.5d0*(xu+xup))
        yv = yq * (yd+0.5d0*(yu+yup))
        zv = zq * (zd+0.5d0*(zu+zup))
        vterm = term * (xq*xq + yq*yq + zq*zq + 2.0d0*(xv+yv+zv)
     &                      + xu*xup + yu*yup + zu*zup
     &                      + xd*(xd+xu+xup) + yd*(yd+yu+yup)
     &                      + zd*(zd+zu+zup))
        vir(1,1) = vir(1,1) + 2.0d0*term*(xq*xq+xv) + vterm
        vir(2,1) = vir(2,1) + 2.0d0*term*(xq*yq+xv)
        vir(3,1) = vir(3,1) + 2.0d0*term*(xq*zq+xv)
        vir(1,2) = vir(1,2) + 2.0d0*term*(yq*xq+yv)
        vir(2,2) = vir(2,2) + 2.0d0*term*(yq*yq+yv) + vterm
        vir(3,2) = vir(3,2) + 2.0d0*term*(yq*zq+yv)
        vir(1,3) = vir(1,3) + 2.0d0*term*(zq*xq+zv)
        vir(2,3) = vir(2,3) + 2.0d0*term*(zq*yq+zv)
        vir(3,3) = vir(3,3) + 2.0d0*term*(zq*zq+zv) + vterm
        if (poltyp .eq. 'DIRECT') then
          vterm = term * (xu*xup+yu*yup+zu*zup)
          vir(1,1) = vir(1,1) + vterm
          vir(2,2) = vir(2,2) + vterm
          vir(3,3) = vir(3,3) + vterm
        end if
      end if
c
c     intermolecular energy is total minus intramolecular part
c
      einter = einter + em - eintra
c      do i = 1, n
c          print*,'dem after Perm x', i, dem(1,i)
c          print*,'dem after Perm y', i, dem(2,i)
c          print*,'dem after Perm z', i, dem(3,i) 
c      end do

      return
      end
