c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine mechanic  --  initialize molecular mechanics  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mechanic" sets up needed parameters for the potential energy
c     calculation and reads in many of the user selectable options
c
c
      subroutine mechanic
      implicit none
      include 'cutoff.i'
      include 'inform.i'
      include 'iounit.i'
      include 'potent.i'
      include 'vdwpot.i'
      include 'ewald.i'
c
c
c     set the bonded connectivity lists and active atoms
c
      call attach
      call active
c
c     find bonds, angles, torsions, bitorsions and small rings
c
      call bonds
      call angles
      call torsions
      call bitors
      call rings
c
c     find unit cell type, lattice parameters and cutoff values
c
      call unitcell
      call lattice
      call polymer
      call cutoffs_serial
c
c     setup needed for potential energy smoothing methods
c
      call flatten
c
c     get the force field parameters and assign atom types
c
      call field
      call katom
c
c     assign atoms to molcules and set the atom groups
c
      call molecule
      call cluster
c
c     find any pisystem atoms, bonds and torsional angles
c
      call orbital
      call cutoffs2_serial
c
c     assign bond, angle and cross term potential parameters
c
      if (use_bond .or. use_strbnd .or. use_strtor
     &    .or. (use_vdw .and. vdwtyp.eq.'MM3-HBOND'))  call kbond
      if (use_angle .or. use_strbnd .or. use_angang)  call kangle
      if (use_strbnd)  call kstrbnd
      if (use_urey)  call kurey
      if (use_angang)  call kangang
c
c     assign out-of-plane deformation potential parameters
c
      if (use_angle .or. use_opbend)  call kopbend
      if (use_angle .or. use_opdist)  call kopdist
      if (use_improp)  call kimprop
      if (use_imptor)  call kimptor
c
c     assign torsion and torsion cross term potential parameters
c
      if (use_tors .or. use_strtor .or. use_tortor)  call ktors
      if (use_pitors)  call kpitors
      if (use_strtor)  call kstrtor
      if (use_tortor)  call ktortor
c
c     assign van der Waals and electrostatic potential parameters
c
      if (use_vdw .or. use_solv)  call kvdw
      if (use_charge .or. use_chgdpl .or. use_solv)  call kcharge
      if (use_dipole .or. use_chgdpl)  call kdipole

c      if (use_mpole .or. use_polar .or.
c     &    use_solv .or. use_rxnfld) call kmpole
      if (use_mpole .or. use_polar .or.
     &    use_solv .or. use_rxnfld) call kmpole_mod

c      if (use_polar .or. use_solv)  call kpolar
      if (use_polar .or. use_solv)  call kpolar_mod
      if (use_ewald)  call kewald
       
c      print*,"aewaldPerm",aewaldPerm
c
c     assign solvation, metal, pisystem and restraint parameters
c
      if (use_solv)  call ksolv
      if (use_metal) call kmetal
      if (use_orbit) call korbit
      if (use_geom)  call kgeom
c
c     set hybrid parameter values for free energy perturbation
c
      call mutate
c
c     quit if essential parameter information is missing
c
      if (abort) then
         write (iout,10)
   10    format (/,' MECHANIC  --  Some Required Potential Energy',
     &              ' Parameters are Undefined')
         call fatal
      end if
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine mechanic  --  initialize molecular mechanics  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mechanic" sets up needed parameters for the potential energy
c     calculation and reads in many of the user selectable options
c
c
      subroutine mechanic_parallel
      implicit none
      include 'cutoff.i'
      include 'inform.i'
      include 'iounit.i'
      include 'potent.i'
      include 'vdwpot.i'
      include 'ewald.i'
c
c
c     set the bonded connectivity lists and active atoms
c
      call attach
      call active
c
c     find bonds, angles, torsions, bitorsions and small rings
c
c      call bonds
c      call angles
c      call torsions
c      call bitors
c      call rings
c
c     find unit cell type, lattice parameters and cutoff values
c
      call unitcell
      call lattice
      call polymer
c      call cutoffs
      call cutoffs_parallel

c
c     setup needed for potential energy smoothing methods
c

      call flatten

c
c     get the force field parameters and assign atom types
c
      call field
      call katom
c
c     assign atoms to molcules and set the atom groups
c

c      call molecule
      call cluster
c
c     find any pisystem atoms, bonds and torsional angles
c
c      call orbital
c      call cutoffs2
      call cutoffs2_parallel
c
c     assign bond, angle and cross term potential parameters
c
c      if (use_bond .or. use_strbnd .or. use_strtor
c     &    .or. (use_vdw .and. vdwtyp.eq.'MM3-HBOND'))  call kbond
c      if (use_angle .or. use_strbnd .or. use_angang)  call kangle
c      if (use_strbnd)  call kstrbnd
c      if (use_urey)  call kurey
c      if (use_angang)  call kangang
c
c     assign out-of-plane deformation potential parameters
c
c      if (use_angle .or. use_opbend)  call kopbend
c      if (use_angle .or. use_opdist)  call kopdist
c      if (use_improp)  call kimprop
c      if (use_imptor)  call kimptor
c
c     assign torsion and torsion cross term potential parameters
c
c      if (use_tors .or. use_strtor .or. use_tortor)  call ktors
c      if (use_pitors)  call kpitors
c      if (use_strtor)  call kstrtor
c      if (use_tortor)  call ktortor
c
c     assign van der Waals and electrostatic potential parameters
c
      if (use_vdw .or. use_solv)  call kvdw
c      if (use_charge .or. use_chgdpl .or. use_solv)  call kcharge
c      if (use_dipole .or. use_chgdpl)  call kdipole

ccc      if (use_mpole .or. use_polar .or.
ccc     &    use_solv .or. use_rxnfld) call kmpole

      if (use_mpole .or. use_polar .or.
     &    use_solv .or. use_rxnfld) call kmpole_mod

ccc      if (use_polar .or. use_solv)  call kpolar
      if (use_polar .or. use_solv)  call kpolar_mod

c      if (use_ewald)  call kewald
      if (use_ewald) call ewaldcof (aewaldPerm,ewaldcut)

      if (use_ewald3b)  call kewald3b

       
c      print*,"aewaldPerm",aewaldPerm
c
c     assign solvation, metal, pisystem and restraint parameters
c
c      if (use_solv)  call ksolv
c      if (use_metal) call kmetal
c      if (use_orbit) call korbit
c      if (use_geom)  call kgeom
c
c     set hybrid parameter values for free energy perturbation
c
c      call mutate
c
c     quit if essential parameter information is missing
c
      if (abort) then
         write (iout,10)
   10    format (/,' MECHANIC  --  Some Required Potential Energy',
     &              ' Parameters are Undefined')
         call fatal
      end if
      return
      end
c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine mechanic  --  initialize molecular mechanics  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mechanic" sets up needed parameters for the potential energy
c     calculation and reads in many of the user selectable options
c
c
      subroutine mechanic_pme
      implicit none
      include 'cutoff.i'
      include 'inform.i'
      include 'iounit.i'
      include 'potent.i'
      include 'vdwpot.i'
      include 'ewald.i'
c
c
c     set the bonded connectivity lists and active atoms
c
      call attach
      call active
c
c     find bonds, angles, torsions, bitorsions and small rings
c
      call bonds
      call angles
      call torsions
      call bitors
      call rings
c
c     find unit cell type, lattice parameters and cutoff values
c
      call unitcell
      call lattice
      call polymer
      call cutoffs_serial
c
c     setup needed for potential energy smoothing methods
c
      call flatten
c
c     get the force field parameters and assign atom types
c
      call field
      call katom
c
c     assign atoms to molcules and set the atom groups
c
      call molecule
      call cluster
c
c     find any pisystem atoms, bonds and torsional angles
c
      call orbital
      call cutoffs2_serial
c
c     assign bond, angle and cross term potential parameters
c
      if (use_bond .or. use_strbnd .or. use_strtor
     &    .or. (use_vdw .and. vdwtyp.eq.'MM3-HBOND'))  call kbond
      if (use_angle .or. use_strbnd .or. use_angang)  call kangle
      if (use_strbnd)  call kstrbnd
      if (use_urey)  call kurey
      if (use_angang)  call kangang
c
c     assign out-of-plane deformation potential parameters
c
      if (use_angle .or. use_opbend)  call kopbend
      if (use_angle .or. use_opdist)  call kopdist
      if (use_improp)  call kimprop
      if (use_imptor)  call kimptor
c
c     assign torsion and torsion cross term potential parameters
c
      if (use_tors .or. use_strtor .or. use_tortor)  call ktors
      if (use_pitors)  call kpitors
      if (use_strtor)  call kstrtor
      if (use_tortor)  call ktortor
c
c     assign van der Waals and electrostatic potential parameters
c
      if (use_vdw .or. use_solv)  call kvdw
      if (use_charge .or. use_chgdpl .or. use_solv)  call kcharge
      if (use_dipole .or. use_chgdpl)  call kdipole

c      if (use_mpole .or. use_polar .or.
c     &    use_solv .or. use_rxnfld) call kmpole
      if (use_mpole .or. use_polar .or.
     &    use_solv .or. use_rxnfld) call kmpole_mod

c      if (use_polar .or. use_solv)  call kpolar
      if (use_polar .or. use_solv)  call kpolar_mod
      if (use_ewald)  call kewald_pme 
       
c      print*,"aewaldPerm",aewaldPerm
c
c     assign solvation, metal, pisystem and restraint parameters
c
      if (use_solv)  call ksolv
      if (use_metal) call kmetal
      if (use_orbit) call korbit
      if (use_geom)  call kgeom
c
c     set hybrid parameter values for free energy perturbation
c
      call mutate
c
c     quit if essential parameter information is missing
c
      if (abort) then
         write (iout,10)
   10    format (/,' MECHANIC  --  Some Required Potential Energy',
     &              ' Parameters are Undefined')
         call fatal
      end if
      return
      end
