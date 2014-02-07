c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  program analyze_parallel_bcast_evenload                ##
c     ##                                                         ##
c     #############################################################
c
c
c     "analyze_parallel" performs several iterations of an 
c     energy/gradient/viral calculation with an MPI-parallelized
c     implementation of the 3-body approximation of the polarization 
c     component of the energy/gradient/virial calculation.  This is  
c     essentially a format for testing the parallelized polarization 
c     energy/gradient/virial code for subsequent implementation in the
c     molecular dynamics routine, dynamic.f
      program analyze_parallel_bcast_evenload_smallsend_nocalcmaster
      implicit none
c      include 'sizes.i'
      include 'sizes3b.i'
c      include 'atoms.i'
      include 'atoms3b.i'
      include 'files.i'
      include 'inform.i'
      include 'iounit.i'
      include 'mpif.h'
      include 'analyz.i'
      include 'bound.i'
      include 'cutoff.i'
      include 'energi.i'
      include 'group3b.i'
c      include 'inter.i'
      include 'potent.i'
      include 'vdwpot.i'
      include 'deriv.i'
      include 'mpole3b.i'
      include 'mpole.i'
      include 'mpole3b_3.i'
      include 'mpole3b_2.i'
      include 'virial.i'
      include 'molcul.i'
      include 'molcul3b.i'
      include 'neigh.i'
      include 'neigh3b.i'
      include 'neigh3b_recv.i'
      logical doanalysis
      integer i1,i3,ii,k,moli2,moli3,count3,iter3
      real*8 energy,step
      real*8 cutoff
      real*8 efindiff(10),analyticg(10)
      real*8 epArr(3)
      integer cntr 
      integer i,j,ixyz,kk,tottriples,ntriples
      integer frame
      integer freeunit
      integer trimtext
      integer list(20)
      real*8 wall,cpu
      logical dosystem,doparam
      logical doenergy,doatom
      logical dolarge,dodetail
      logical doprops,doconect
      logical exist
      logical, allocatable :: active(:)
      character*1 letter
      character*120 record
      character*120 string
      character*120 xyzfile
      integer ierr,taskid,numtasks,master
      integer l1,i2
      integer k1,k2
      integer moli1
      real*8 ep3bt,virep3bt(3,3),ep3bt_tot2,virep3bt_tot2(3,3)
      real*8 ep3bt_tot3,virep3bt_tot3(3,3)
      real*8 ep3b2,ep3b3,virep3b2(3,3),virep3b3(3,3)
      real*8, allocatable :: dep3bt(:,:)
      real*8, allocatable :: dep3bt_tot2(:,:)
      real*8, allocatable :: dep3b2(:,:)
      real*8, allocatable :: dep3bt_tot3(:,:)
      real*8, allocatable :: dep3b3(:,:)
      integer offset,remainder,start
      integer localsum,totalsum,localsum_tot
      integer stat(MPI_STATUS_SIZE)
c      integer,allocatable :: mol3new(:)
      character*6 mode
      integer s1,s2,s3,s4,totsize1,totsize2
      integer, allocatable :: buf1(:)
      integer, allocatable :: buf2(:)


      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,taskid,ierr)
      call mpi_comm_size(mpi_comm_world,numtasks,ierr)

c
c    Set 0th taskid to master task
c
      master=0


      if(taskid.eq.master) then
      call settime
      call promo
      end if
c      print*,"taskid=",taskid

c     
c     All tasks perform setup of 
c     of the structure and mechanics calculation
c

      call initial
      call getxyz
c      call mechanic
      if(taskid.eq.master) then
         call mechanic
      end if
c       call mpi_barrier(mpi_comm_world,ierr)

         call mpi_bcast(nmol,1,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(imol,2*maxatm,mpi_integer,master,
     &   mpi_comm_world,ierr)
c         print*,"After molecule bcast"

      if(taskid.ne.master) then
         call mechanic_parallel
      end if

c          call mpi_bcast(ipole,maxatm,mpi_integer,master,
c     &   mpi_comm_world,ierr)
c          call mpi_bcast(xaxis,maxatm,mpi_integer,master,
c     &   mpi_comm_world,ierr)
c          call mpi_bcast(yaxis,maxatm,mpi_integer,master,
c     &   mpi_comm_world,ierr)
c          call mpi_bcast(zaxis,maxatm,mpi_integer,master,
c     &   mpi_comm_world,ierr)
c          call mpi_bcast(npole,1,mpi_integer,master,
c     &   mpi_comm_world,ierr)
c          call mpi_bcast(polaxe,8*maxatm,mpi_character,master,
c     &   mpi_comm_world,ierr)
c         print*,"After mpole var bcast"
c     All tasks allocate 2- and 3-body gradient of polarization
c     energy(dep3b2 and dep3b3, as well as the temporary running 
c     totals of these, dep3bt_tot2 and dep3bt_tot3.  

c      allocate (dep3bt(3,npole))
      if(taskid.eq.master) then
        allocate (dep3b2(3,npole))
        allocate (dep3b3(3,npole))
      end if
      if(taskid.ne.master) then
      allocate (dep3bt_tot2(3,npole))
      allocate (dep3bt_tot3(3,npole))
      allocate (dep3bt(3,npole))
      end if
c     An array storing the nonzero elements of the 3-body neighbor list.
c      allocate (mol3new(nmol))
c
c     get the desired types of analysis to be performed
c
      call nextarg (string,exist)
      if (.not. exist) then
         write (iout,10)
   10    format (/,' The TINKER Analysis Facility can Provide :',
     &           //,' General System and Force Field Information [G]',
     &           /,' Force Field Parameters for Interactions [P]',
     &           /,' Total Potential Energy and its Components [E]',
     &           /,' Energy Breakdown over Each of the Atoms [A]',
     &           /,' List of the Large Individual Interactions [L]',
     &           /,' Details for All Individual Interactions [D]',
     &           /,' Electrostatic, Inertial & Virial Properties [M]',
     &           /,' Connectivity Lists for Each of the Atoms [C]')
   20    continue
         write (iout,30)
   30    format (/,' Enter the Desired Analysis Types',
     &              ' [G,P,E,A,L,D,M] :  ',$)
         read (input,40,err=20)  string
   40    format (a120)
      end if
c
c     set option control flags based desired analysis types
c
      dosystem = .false.
      doparam = .false.
      doenergy = .false.
      doatom = .false.
      dolarge = .false.
      dodetail = .false.
      doprops = .false.
      doconect = .false.
      call upcase (string)
      do i = 1, trimtext(string)
         letter = string(i:i)
         if (letter .eq. 'G')  dosystem = .true.
         if (letter .eq. 'P')  doparam = .true.
         if (letter .eq. 'E')  doenergy = .true.
         if (letter .eq. 'A')  doatom = .true.
         if (letter .eq. 'L')  dolarge = .true.
         if (letter .eq. 'D')  dodetail = .true.
         if (letter .eq. 'M')  doprops = .true.
         if (letter .eq. 'C')  doconect = .true.
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (active(n))
c
c     get the list of atoms for which output is desired
c
      if (doatom .or. doparam) then
         do i = 1, 20
            list(i) = 0
         end do
         if (exist) then
            do i = 1, 20
               call nextarg (string,exist)
               if (.not. exist)  goto 50
               read (string,*,err=50,end=50)  list(i)
            end do
   50       continue
         else
            write (iout,60)
   60       format (/,' List Atoms for which Output is Desired',
     &                 ' [ALL] :  '/,'    >  ',$)
            read (input,70)  record
   70       format (a120)
            read (record,*,err=80,end=80)  (list(i),i=1,20)
   80       continue
         end if
         do i = 1, n
            active(i) = .true.
         end do
         i = 1
         do while (list(i) .ne. 0)
            if (i .eq. 1) then
               do j = 1, n
                  active(j) = .false.
               end do
            end if
            if (list(i) .gt. 0) then
               active(list(i)) = .true.
               i = i + 1
            else
               do j = abs(list(i)), abs(list(i+1))
                  active(j) = .true.
               end do
               i = i + 2
            end if
         end do
      end if
c
c     setup to write out the large individual energy terms
c
      if (dolarge) then
         verbose = .true.
      end if
c
c     setup to write out all of the individual energy terms
c
      if (dodetail) then
         doenergy = .true.
         debug = .true.
         verbose = .true.
      else
         debug = .false.
      end if
c
c     reopen the coordinates file and read the first structure
c
      frame = 0
      ixyz = freeunit ()
      xyzfile = filename
      call suffix (xyzfile,'xyz','old')
      open (unit=ixyz,file=xyzfile,status ='old')
      rewind (unit=ixyz)
      call readxyz (ixyz)
c
c     get info on the molecular system and force field
c
      if (dosystem)  call systyze
c
c     get parameters used for molecular mechanics potentials
c
      if (doparam .and. doconect) then
         call amberyze (active)
      else if (doparam) then
         call paramyze (active)
      end if


c
c     make the call to compute the potential energy,gradient and virial
c
         if (doenergy .or. doatom .or. dolarge) then
c            call enrgyze()
            doanalysis = .true.
         end if
      step=0.0000001d0

c  Build 2- and 3-body neighbor lists
             offset=int(nmol/(numtasks-1))
             remainder=mod(nmol,numtasks-1)
c      print*,"Before pack size1"
c      call mpi_pack_size(800*offset,mpi_integer,mpi_comm_world,
c     &     s1,ierr)

c      print*,"Before pack size2"
c      call mpi_pack_size(offset,mpi_integer,mpi_comm_world,
c     &     s2,ierr)

c      totsize1=numtasks*(s1+s2+2*mpi_bsend_overhead)
c      allocate (buf1(totsize1))

c      print*,"Before 1st buffer attach"
c      call mpi_buffer_attach(buf1,totsize1,ierr)
c      print*,"After 1st buffer attach"

c      if (taskid.le.remainder-1) then
c       call mpi_pack_size(800,mpi_integer,mpi_comm_world,
c     &     s3,ierr)
c       call mpi_pack_size(1,mpi_integer,mpi_comm_world,
c     &     s4,ierr)
c       totsize2=remainder*(s3+s4+2*mpi_bsend_overhead)
c       allocate (buf2(totsize2))
c       call mpi_buffer_attach(buf2,totsize2,ierr)
c      end if

      if(taskid.eq.master) then
        if(domollst2bod) call molfull2body

        if(domollst3bod) call molfull3bodycobar
c           do i=1,numtasks
c           call mpi_bsend(mollst3mod(1:800,((i-1)*offset+1):i*offset),
c     &      800*offset,mpi_integer,i-1,2*(i-1),mpi_comm_world,ierr)
c           call mpi_bsend(nmollst3mod(((i-1)*offset+1):i*offset),offset,
c     &      mpi_integer,i-1,2*(i-1)+1,mpi_comm_world,ierr)
c           end do
c           print*,"After 1st send"
c           do i=1,remainder
c             call mpi_bsend(mollst3mod(1:800,
c     &      (numtasks*offset+i):(numtasks*offset+i)),800,mpi_integer,
c     &      i-1,2*numtasks+2*(i-1),mpi_comm_world,ierr)
c             call mpi_bsend(nmollst3mod((numtasks*offset+i):
c     &        (numtasks*offset+i)),1,mpi_integer,i-1,
c     &      2*numtasks+2*(i-1)+1,mpi_comm_world,ierr)
c           end do
c          print*,"After 2nd send"
      end if

c      call mpi_pack_size(800,mpi_integer,mpi_comm_world,
c     &     s1,ierr)

c      call mpi_pack_size(1,mpi_integer,mpi_comm_world,
c     &     s2,ierr)

c      totsize1=numtasks*offset*(s1+s2+2*mpi_bsend_overhead)
c      allocate (buf1(totsize1))
c      print*,"totsize1",totsize1

c      if (taskid.le.remainder-1) then
c       call mpi_pack_size(800,mpi_integer,mpi_comm_world,
c     &     s3,ierr)
c       call mpi_pack_size(1,mpi_integer,mpi_comm_world,
c     &     s4,ierr)
c       totsize2=remainder*(s3+s4+2*mpi_bsend_overhead)
c       allocate (buf2(totsize2))
c       print*,"totsize2",totsize2
c       call mpi_buffer_attach(buf2,totsize2,ierr)
c      end if
c       allocate (mollst3mod_recv(800,offset))
c       allocate (nmollst3mod_recv(offset))
       allocate (mollst3mod_recv_small(800,1))
       allocate (nmollst3mod_recv_small(1))
c          call mpi_recv(mollst3mod_recv,800*offset,mpi_integer,master,
c     &      2*taskid,mpi_comm_world,stat,ierr)
c          call mpi_recv(nmollst3mod_recv,offset,mpi_integer,master,
c     &      2*taskid+1,mpi_comm_world,stat,ierr)

c         if(taskid.le.remainder-1) then
c          call mpi_recv(mollst3mod_recv_small,800,mpi_integer,master,
c     &      2*numtasks+2*(taskid),mpi_comm_world,stat,ierr)
c          call mpi_recv(nmollst3mod_recv_small,1,mpi_integer,master,
c     &      2*numtasks+2*(taskid)+1,mpi_comm_world,stat,ierr)
c         end if

c
c   Broadcast neighbor lists two all tasks
c
         call mpi_bcast(nmollst,nmol,mpi_integer,master,
     &   mpi_comm_world,ierr)
c         print*,"After 1st mpi_bcast"
          call mpi_bcast(mollst,100*nmol,mpi_integer,master,
     &   mpi_comm_world,ierr)
c         print*,"After 2nd mpi_bcast"
c         call mpi_bcast(nmollst3mod,nmol,mpi_integer,master,
c     &   mpi_comm_world,ierr)
c         print*,"After 3rd mpi_bcast"
c         call mpi_bcast(mollst3mod,800*nmol,mpi_integer,master,
c     &   mpi_comm_world,ierr)
c         call mpi_barrier(mpi_comm_world,ierr)
c         print*,"After all mpi_bcast"
c
c      Determine elements of 3-body neighbor list that are nonzero
c

c        count3=0
c         do k=1,nmol
c          if(nmollst3mod(k).ne.0) then
c            count3=count3+1
c            mol3new(count3)=k
c          end if
c         end do

      do i3=1,10
        x(3)=x(3)+step
c     Initialize all energies (scalars whose variable names begin w/
c     “e”) and initialize all  gradients (“2-dimensional arrays with
c     variable names beginning  w/ “de”), and virial tensor to 0.0d0
      if(taskid.eq.master) then

         if(doanalysis) then

      eb = 0.0d0
      ea = 0.0d0
      eba = 0.0d0
      eub = 0.0d0
      eaa = 0.0d0
      eopb = 0.0d0
      eopd = 0.0d0
      eid = 0.0d0
      eit = 0.0d0
      et = 0.0d0
      ept = 0.0d0
      ebt = 0.0d0
      ett = 0.0d0
      ev = 0.0d0
      ec = 0.0d0
      ecd = 0.0d0
      ed = 0.0d0
      em = 0.0d0
      ep = 0.0d0
      er = 0.0d0
      es = 0.0d0
      elf = 0.0d0
      eg = 0.0d0
      ex = 0.0d0
      totalsum=0
      ep3b2=0.0d0
      ep3b3=0.0d0
        do i=1,3
          do j=1,3
            virep3b2(i,j)=0.0d0
            virep3b3(i,j)=0.0d0
            vir(i,j)=0.0d0
          end do
        end do

c
c     zero out energy partitioning components for each atom
c
      do i = 1, n
         do j=1,3
            deb(j,i) = 0.0d0
            dea(j,i) = 0.0d0
            deba(j,i) = 0.0d0
            deub(j,i) = 0.0d0
            deaa(j,i) = 0.0d0
            deopb(j,i) = 0.0d0
            deopd(j,i) = 0.0d0
            deid(j,i) = 0.0d0
            deit(j,i) = 0.0d0
            det(j,i) = 0.0d0
            dept(j,i) = 0.0d0
            debt(j,i) = 0.0d0
            dett(j,i) = 0.0d0
            dev(j,i) = 0.0d0
            dec(j,i) = 0.0d0
            decd(j,i) = 0.0d0
            ded(j,i) = 0.0d0
            dem(j,i) = 0.0d0
            dep(j,i) = 0.0d0
            der(j,i) = 0.0d0
            des(j,i) = 0.0d0
            delf(j,i) = 0.0d0
            deg(j,i) = 0.0d0
            dex(j,i) = 0.0d0
            dep3b2(j,i) = 0.0d0
            dep3b3(j,i) = 0.0d0
         end do
         aeb(i) = 0.0d0
         aea(i) = 0.0d0
         aeba(i) = 0.0d0
         aeub(i) = 0.0d0
         aeaa(i) = 0.0d0
         aeopb(i) = 0.0d0
         aeopd(i) = 0.0d0
         aeid(i) = 0.0d0
         aeit(i) = 0.0d0
         aet(i) = 0.0d0
         aept(i) = 0.0d0
         aebt(i) = 0.0d0
         aett(i) = 0.0d0
         aev(i) = 0.0d0
         aec(i) = 0.0d0
         aecd(i) = 0.0d0
         aed(i) = 0.0d0
         aem(i) = 0.0d0
         aep(i) = 0.0d0
         aer(i) = 0.0d0
         aes(i) = 0.0d0
         aelf(i) = 0.0d0
         aeg(i) = 0.0d0
         aex(i) = 0.0d0
      end do
c
c     zero out the total intermolecular energy
c
c      einter = 0.0d0
         end if 
      end if
c
c     maintain any periodic boundary conditions
c
         if(doanalysis) then

      if (use_bounds .and. .not.use_group)  call bounds
c
c     remove any previous use of the replicates method
c
      cutoff = 0.0d0
      call replica (cutoff)
         end if 
c
c     update the pairwise interaction neighbor lists
c
      if(taskid.eq.master) then
         if(doanalysis) then
      if (use_list)  call nblist
c
c     many implicit solvation models require Born radii
c
      if (use_born)  call born
c
c     alter bond and torsion constants for pisystem
c
      if (use_orbit)  call picalc
c
c     call the local geometry energy component routines
c
      if (use_bond)  call ebond1
      if (use_angle)  call eangle1
      if (use_strbnd)  call estrbnd1
      if (use_urey)  call eurey1
      if (use_angang)  call eangang1
      if (use_opbend)  call eopbend1
      if (use_opdist)  call eopdist1
      if (use_improp)  call eimprop1
      if (use_imptor)  call eimptor1
      if (use_tors)  call etors1
      if (use_pitors)  call epitors1
      if (use_strtor)  call estrtor1
      if (use_tortor)  call etortor1
c
c     call the van der Waals energy component routines
c
      if (use_vdw) then
         if (vdwtyp .eq. 'LENNARD-JONES')  call elj1
         if (vdwtyp .eq. 'BUCKINGHAM')  call ebuck1
         if (vdwtyp .eq. 'MM3-HBOND')  call emm3hb1
         if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal1
         if (vdwtyp .eq. 'GAUSSIAN')  call egauss1
      end if
      print*,"After vdw"
c
c     call the electrostatic energy component routines
c
      if (use_charge)  call echarge1
      if (use_chgdpl)  call echgdpl1
      if (use_dipole)  call edipole1
c      if (use_mpole .or. use_polar)  call empole1c_3b
c      if (use_mpole .or. use_polar) then
c             do_empole1c_3b=.true.
c      end if

c    Call permanent electrostatics calculation
        print*,"Before Perm Elec"
        call empole1c_3b_Perm
        print*,"Perm Elec Before Parallel",em

      if (use_rxnfld)  call erxnfld1
c
c     call any miscellaneous energy component routines
c
      if (use_solv)  call esolv1
      if (use_metal)  call emetal1
      if (use_geom)  call egeom1
      if (use_extra)  call extra1

c       tottriples=0
c       if (i3.eq.1) then
c         do k=1,nmol
c          ntriples=nmollst3mod(k)/2
c          print*,"Nmollst3mod",k,nmollst3mod(k)
c          tottriples=tottriples+ntriples
c         end do
c         print*,"Triplecount=",tottriples
c       end if

         end if

      end if

 

c    Broadcast multipoles (variable name "rpole") to all 
c    slave tasks after function called within empole1c_3b_Perm
c    transformes them from the locally defined coordinate frame to 
c    that of the whole system.

c      if(i3.eq.1) then   
         call mpi_bcast(rpole,13*maxatm,mpi_real8,master,
     &   mpi_comm_world,ierr)
c         call mpi_barrier(mpi_comm_world,ierr)
c      end if
c       print*,"After multipole bcast"
c    Begin polarization energy,gradient, and virial calculation.  
c    First, call switch to set up the calculation for the infinite cutoff specified
c    by the 'MPOLE' flag.  Call barrier to ensure that all tasks have
c    completed this subroutine call before proceeding.

      mode = 'MPOLE'
      call switch (mode)
c         call mpi_barrier(mpi_comm_world,ierr)

         if(doanalysis) then
c     Ininitialize running totals for 2- and 3-body polarization
c     energy,gradient, and virial
c
               if(taskid.ne.master) then
                ep3bt_tot2=0.0d0
                ep3bt_tot3=0.0d0
                do i = 1, npole
                   do j = 1, 3
                      dep3bt_tot2(j,i) = 0.0d0
                      dep3bt_tot3(j,i) = 0.0d0
                   end do
                end do

                do i=1,3
                   do j=1,3
                      virep3bt_tot2(i,j)=0.0d0
                      virep3bt_tot3(i,j)=0.0d0
                   end do
                end do
               end if
c    Determine partitioning of 2-body polarization energy calculations
c    among tasks.  Each task gets at least 'offset' number of tasks,
c    where nmol is the total number of molecules.  If taskid is less
c    than the 'remainder', that particular task will do one more iteration,
c    in addition to the number of iterations specified by 'offset'.
c    This is for load balancing.

c               call mpi_buffer_attach(buf1,totsize1,ierr)
c               print*,"After buf1 attach"
c               if (taskid.le.remainder-1) then
c                 call mpi_buffer_attach(buf2,totsize2,ierr)
c               end if
c               print*,"After buf2 attach"

               if(taskid.eq.master) then
                 do i=1,numtasks-1
                   do j=1,offset
                   call mpi_send(mollst3mod(1:800,(i-1)*offset+j:
     &           (i-1)*offset+j),800,mpi_integer,i,
     &           2*offset*(i-1)+2*(j-1),mpi_comm_world,ierr)
                   call mpi_send(nmollst3mod( (i-1)*offset+j:
     &           (i-1)*offset+j),1,mpi_integer,i,
     &           2*offset*(i-1)+2*(j-1)+1,mpi_comm_world,ierr)
                   end do
                 end do
                 do i=1,remainder
                 call mpi_send(mollst3mod(1:800,((numtasks-1)*offset+i):
     &            ((numtasks-1)*offset+i)),800,mpi_integer,i,
     &             2*(numtasks-1)*offset+2*(i-1),mpi_comm_world,ierr)
                 call mpi_send(nmollst3mod(((numtasks-1)*offset+i):
     &            ((numtasks-1)*offset+i)),1,mpi_integer,i,
     &              2*(numtasks-1)*offset+2*(i-1)+1,mpi_comm_world,ierr)
                 end do
               end if

c             offset=int(nmol/numtasks)
c             remainder=mod(nmol,numtasks)

c      Each task then goes through its number of iterations (determined
c      above).
c      Within each iteration, Innerloop2 is called (needs the body index
c      (mol1) as input), returning ep3bt, virep3bt, dep3bt which are
c      summed within that task. Once the totals for each task are calculated
c      they are mpi_reduced.


             if((taskid.le.remainder).and.(taskid.ne.master)) then
               start=(taskid-1)*offset+1
               do moli1 =start,start+offset-1
                  call Innerloop2(moli1,ep3bt,virep3bt,dep3bt)
                 do i=1,3
                   do j=1,3
                 virep3bt_tot2(i,j)=virep3bt_tot2(i,j)+virep3bt(i,j)
                   end do
                 end do

                 do i=1,npole
                    do j=1,3
                 dep3bt_tot2(j,i)=dep3bt_tot2(j,i)+dep3bt(j,i)
                    end do
                 end do

                 ep3bt_tot2=ep3bt_tot2+ep3bt                 
               end do 

                moli1=(numtasks-1)*offset+taskid
                  call Innerloop2(moli1,ep3bt,virep3bt,dep3bt)
                 do i=1,3
                   do j=1,3
                   virep3bt_tot2(i,j)=virep3bt_tot2(i,j)+virep3bt(i,j)
                   end do
                 end do
                 do i=1,npole
                    do j=1,3
                    dep3bt_tot2(j,i)=dep3bt_tot2(j,i)+dep3bt(j,i)
                    end do
                 end do
                ep3bt_tot2=ep3bt_tot2+ep3bt
               print*,"ep3bt_tot2",ep3bt_tot2

                  call mpi_reduce(ep3bt_tot2,ep3b2,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt_tot2,dep3b2,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt_tot2,virep3b2,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
             end if
c             else if(taskid.gt.remainder) then
             if(taskid.gt.remainder) then
               start=(taskid-1)*offset+1
               do moli1 =start,start+offset-1
                  call Innerloop2(moli1,ep3bt,virep3bt,dep3bt)

                 do i=1,3
                   do j=1,3
                    virep3bt_tot2(i,j)=virep3bt_tot2(i,j)+virep3bt(i,j)
                   end do
                 end do

                 do i=1,npole
                    do j=1,3
                    dep3bt_tot2(j,i)=dep3bt_tot2(j,i)+dep3bt(j,i)
                    end do
                 end do

                ep3bt_tot2=ep3bt_tot2+ep3bt
               end do 
               print*,"ep3bt_tot2",ep3bt_tot2
                  call mpi_reduce(ep3bt_tot2,ep3b2,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt_tot2,dep3b2,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt_tot2,virep3b2,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
             end if
c    Determine partitioning of 3-body polarization energy calculations
c    among tasks.  Each task gets at least 'offset' number of tasks,
c    where count3 is the total number of molecules whose 3-body neighbor
c    list count is nonzero.  If taskid is less
c    than the 'remainder', that particular task will do one more
c    iteration, in addition to the number of iterations specified by 'offset'.
c    This is for load balancing.

c             offset=int(count3/numtasks)
c             remainder=mod(count3,numtasks)
c      Each task then goes through its number of iterations (determined
c      above).
c      Within each iteration, Innerloop3 is called (needs the body index
c      (mol1) as input), returning ep3bt, virep3bt, dep3bt which are
c      summed within that task. Once the totals for each task are
c      calculated they are mpi_reduced.

c             if(taskid.le.remainder) then
             if((taskid.le.remainder).and.(taskid.ne.master)) then
               start=(taskid-1)*offset+1
               iter3=0
               cntr=0
               do moli1 =start,start+offset-1
c               do k1 =start,start+offset-1
c                 moli1=mol3new(k1)
                cntr=cntr+1
c  2*offset*(i-1)+2*(j-1)+1
                call mpi_recv(mollst3mod_recv_small,800,mpi_integer,
     &          master,2*offset*(taskid-1)+2*(cntr-1),mpi_comm_world,
     &          stat,ierr)
                call mpi_recv(nmollst3mod_recv_small,1,mpi_integer,
     &          master,2*offset*(taskid-1)+2*(cntr-1)+1,mpi_comm_world,
     &          stat,ierr)

c                 iter3=iter3+1
                iter3=1                
c                 call Innerloop3(iter3,moli1,ep3bt,virep3bt,dep3bt)
                call Innerloop3_small(iter3,moli1,ep3bt,virep3bt,dep3bt)

                 do i=1,3
                   do j=1,3
                 virep3bt_tot3(i,j)=virep3bt_tot3(i,j)+virep3bt(i,j)
                   end do
                 end do

                 do i=1,npole
                    do j=1,3
                 dep3bt_tot3(j,i)=dep3bt_tot3(j,i)+dep3bt(j,i)
                    end do
                 end do

                 ep3bt_tot3=ep3bt_tot3+ep3bt
               end do

                iter3=1
c                moli1=numtasks*offset+taskid+1
cc2*(numtasks-1)*offset+2*(i-1)+1
                moli1=(numtasks-1)*offset+taskid
cccc                k1=numtasks*offset+taskid+1
cccc                moli1=mol3new(k1)
                call mpi_recv(mollst3mod_recv_small,800,mpi_integer,
     &          master,2*(numtasks-1)*offset+2*(taskid-1),
     &          mpi_comm_world,stat,ierr)
                call mpi_recv(nmollst3mod_recv_small,1,mpi_integer,
     &          master,2*(numtasks-1)*offset+2*(taskid-1)+1,
     &          mpi_comm_world,stat,ierr)

cccc 2*numtasks*offset+2*(i-1)

                call Innerloop3_small(iter3,moli1,ep3bt,virep3bt,dep3bt)
                 do i=1,3
                   do j=1,3
                   virep3bt_tot3(i,j)=virep3bt_tot3(i,j)+virep3bt(i,j)
                   end do
                 end do
                 do i=1,npole
                    do j=1,3
                    dep3bt_tot3(j,i)=dep3bt_tot3(j,i)+dep3bt(j,i)
                    end do
                 end do
                ep3bt_tot3=ep3bt_tot3+ep3bt
                print*,"ep3bt_tot3 first part",ep3bt_tot3
                  call mpi_reduce(ep3bt_tot3,ep3b3,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt_tot3,dep3b3,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt_tot3,virep3b3,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
             end if
             if(taskid.gt.remainder) then 
               start=(taskid-1)*offset+1
               iter3=0
               cntr=0
               do moli1 =start,start+offset-1
c               do k1=start,start+offset-1
c                  moli1=mol3new(k1)
                cntr=cntr+1
                call mpi_recv(mollst3mod_recv_small,800,mpi_integer,
     &          master,2*offset*(taskid-1)+2*(cntr-1),mpi_comm_world,
     &          stat,ierr)
                call mpi_recv(nmollst3mod_recv_small,1,mpi_integer,
     &          master,2*offset*(taskid-1)+2*(cntr-1)+1,mpi_comm_world,
     &          stat,ierr)

c                  iter3=iter3+1
                 iter3=1
c                  call Innerloop3(iter3,moli1,ep3bt,virep3bt,dep3bt)
                call Innerloop3_small(iter3,moli1,ep3bt,virep3bt,dep3bt)


                 do i=1,3
                   do j=1,3
                    virep3bt_tot3(i,j)=virep3bt_tot3(i,j)+virep3bt(i,j)
                   end do
                 end do

                 do i=1,npole
                    do j=1,3
                    dep3bt_tot3(j,i)=dep3bt_tot3(j,i)+dep3bt(j,i)
                    end do
                 end do

                ep3bt_tot3=ep3bt_tot3+ep3bt
               end do
                print*,"ep3bt_tot3 second part",ep3bt_tot3

                  call mpi_reduce(ep3bt_tot3,ep3b3,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt_tot3,dep3b3,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt_tot3,virep3b3,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
             end if

c     Call barrier to ensure that all tasks have completed their portion
c     of the polarization energy calculation.     
c
c        call mpi_buffer_detach(buf1,totsize1,ierr)
c        call mpi_buffer_detach(buf2,totsize2,ierr)
           call mpi_barrier(mpi_comm_world,ierr)          

c     Master task sums 2- and 3-body polarization energy, gradient, and
c     virial, and sums total gradient, energy, and virial due to the
c    other noncovalent components, as well as the covalent energies.
           if(taskid.eq.master) then
             ep = ep+ep3b2+ep3b3
             print*,"In master ep3b2",ep3b2
             print*,"In master ep3b3",ep3b3

             do i = 1, npole
               do j = 1, 3
                 dep(j,i) =dep(j,i)+dep3b2(j,i)+dep3b3(j,i)
               end do
             end do

             do i=1,3
               do j=1,3
                 vir(i,j)=vir(i,j)+virep3b2(i,j)+virep3b3(i,j)
               end do
             end do

      esum = eb + ea + eba + eub + eaa + eopb + eopd + eid + eit
     &          + et + ept + ebt + ett + ev + ec + ecd + ed + em
     &          + ep + er + es + elf + eg + ex
      energy = esum
      do i = 1, n
         do j = 1, 3
            desum(j,i) = deb(j,i) + dea(j,i) + deba(j,i)
     &                      + deub(j,i) + deaa(j,i) + deopb(j,i)
     &                      + deopd(j,i) + deid(j,i) + deit(j,i)
     &                      + det(j,i) + dept(j,i) + debt(j,i)
     &                      + dett(j,i) + dev(j,i) + dec(j,i)
     &                      + decd(j,i) + ded(j,i) + dem(j,i)
     &                      + dep(j,i) + der(j,i) + des(j,i)
     &                      + delf(j,i) + deg(j,i) + dex(j,i)
         end do
      end do

      efindiff(i3)=esum
      analyticg(i3)=desum(1,3)
      if(i3.eq.1) then
        do i=1,n
         print*,"Dep i x",dep(1,i)
         print*,"Dep i y",dep(2,i)
         print*,"Dep i z",dep(3,i)
        end do
      end if
      print*,"Perm Elec Eng, Polarization Eng",em,ep
      print*,"Total Potential Energy :",energy,"Kcal/mole"

           end if  

         end if
      end do


c
c     perform deallocation of some local arrays
c
      deallocate (active)
c
c     perform any final tasks before program exit
c
      close (unit=ixyz)
      if (dodetail)  debug = .false.
      call final

      if(taskid.eq.master) then
      call gettime (wall,cpu)
      print*,"Finite Diff, AnalyticGrad",
     & (efindiff(3)-efindiff(1))/2.0d0/step,analyticg(2)
      print*,"Wall Time=", wall, "CPU Time=",cpu
      end if

c      deallocate(dep3bt)
c      deallocate(dep3bt_tot2)
c      deallocate(dep3bt_tot3)
      if(taskid.ne.master) then
        deallocate(dep3bt)
        deallocate(dep3bt_tot2)
        deallocate(dep3bt_tot3)
      end if
      if(taskid.eq.master) then
        deallocate(dep3b2)
        deallocate(dep3b3)
      end if
c      if (taskid.le.remainder-1) then
c         deallocate(buf2)
c      end if
c       deallocate(buf1)
c       deallocate(mollst3mod_recv)
c       deallocate(nmollst3mod_recv)
       deallocate(mollst3mod_recv_small)
       deallocate(nmollst3mod_recv_small)
c      deallocate(mol3new)
      call mpi_finalize(ierr)
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine systyze  --  system & force field information  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "systyze" is an auxiliary routine for the analyze program
c     that prints general information about the molecular system
c     and the force field model
c
c
      subroutine systyze
      implicit none
      include 'sizes3b.i'
      include 'sizes.i'
      include 'atoms.i'
      include 'atoms3b.i'
      include 'bound.i'
      include 'boxes.i'
      include 'cutoff.i'
      include 'ewald.i'
      include 'fields.i'
      include 'iounit.i'
      include 'molcul.i'
      include 'molcul3b.i'
      include 'pme.i'
      include 'potent.i'
      include 'units.i'
      include 'vdwpot.i'
      integer i
      real*8 dens
      character*20 value
      character*20 label(5)
c
c
c     info on number of atoms, molecules and mass
c
      if (n .ne. 0) then
         write (iout,10)  n,nmol,totmass
   10    format (/,' Overall System Contents :',
     &           //,' Number of Atoms',25x,i8,
     &           /,' Number of Molecules',21x,i8,
     &           /,' Total System Mass',19x,f12.4)
         if (use_bounds) then
            dens = (1.0d24/volbox) * (totmass/avogadro)
            write (iout,20)  dens
   20       format (' System Density',22x,f12.4)
         end if
      end if
c
c     periodic box dimensions and crystal lattice type
c
      if (use_bounds) then
         value = 'ORTHOGONAL'
         if (monoclinic)  value = 'MONOCLINIC'
         if (triclinic)  value = 'TRICLINIC'
         if (octahedron)  value = 'TRUNCATED OCTAHEDRON'
         call justify (value)
         write (iout,30)  xbox,ybox,zbox,alpha,beta,gamma,volbox,value
   30    format (/,' Periodic Boundary Conditions :',
     &           //,' a-Axis Length',23x,f12.4,
     &           /,' b-Axis Length',23x,f12.4,
     &           /,' c-Axis Length',23x,f12.4,
     &           /,' Alpha Angle',25x,f12.4,
     &           /,' Beta Angle',26x,f12.4,
     &           /,' Gamma Angle',25x,f12.4,
     &           /,' Cell Volume',25x,f12.4,
     &           /,' Lattice Type',16x,a20)
         if (spacegrp .ne. '          ') then
            value = spacegrp
            call justify (value)
            write (iout,40)  value
   40       format (' Space Group',17x,a20)
         end if
      end if
c
c     info on force field potential energy function
c
      value = forcefield
      call justify (value)
      write (iout,50)  value
   50 format (/,' Force Field Name :',10x,a20)
c
c     details of vdw potential energy functional form
c
      if (use_vdw) then
         label(1) = vdwtyp
         label(2) = radtyp
         label(3) = radsiz
         label(4) = radrule
         label(5) = epsrule
         do i = 1, 5
            call justify (label(i))
         end do
         write (iout,60)  (label(i),i=1,5)
   60    format (/,' VDW Function',16x,a20,
     &           /,' Size Descriptor',13x,a20,
     &           /,' Size Unit Type',14x,a20,
     &           /,' Size Combining Rule',9x,a20,
     &           /,' Well Depth Rule',13x,a20)
         if (vdwcut .le. 1000.0d0) then
            write (iout,70)  vdwcut
   70       format (' VDW Cutoff',26x,f12.4)
         end if
      end if
c
c     details of electrostatic energy functional form
c
      if (use_charge .or. use_dipole .or. use_mpole .or. use_polar) then
         write (iout,80)
   80    format ()
      end if
      if (use_charge) then
         value = 'PARTIAL CHARGE'
         call justify (value)
         write (iout,90)  value
   90    format (' Electrostatics',14x,a20)
      end if
      if (use_dipole) then
         value = 'BOND DIPOLE'
         call justify (value)
         write (iout,100)  value
  100    format (' Electrostatics',14x,a20)
      end if
      if (use_mpole) then
         value = 'ATOMIC MULTIPOLE'
         call justify (value)
         write (iout,110)  value
  110    format (' Electrostatics',14x,a20)
      end if
      if (use_polar) then
         value = 'INDUCED DIPOLE'
         call justify (value)
         write (iout,120)  value
  120    format (' Electrostatics',14x,a20)
      end if
c
c     details of particle mesh Ewald calculation
c
c      if (use_ewald) then
c         value = boundary
c         call justify (value)
c         write (iout,130)  aewald,ewaldcut,nfft1,nfft2,nfft3,
c     &                    bsorder,value
c  130    format (/,' Particle Mesh Ewald :',
c     &           //,' Ewald Coefficient',19x,f12.4,
c     &           /,' Real-Space Cutoff',19x,f12.4,
c     &           /,' Grid Dimensions',21x,3i4,
c     &           /,' B-Spline Order',26x,i8,
c     &           /,' Boundary Condition',10x,a20)
c      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine paramyze  --  force field parameter analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "paramyze" prints the force field parameters used in the
c     computation of each of the potential energy terms
c
c
      subroutine paramyze (active)
      implicit none
      include 'sizes3b.i'
      include 'sizes.i'
      include 'angang.i'
      include 'angle.i'
      include 'angpot.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'atoms3b.i'
      include 'bitor.i'
      include 'bond.i'
      include 'charge.i'
      include 'dipole.i'
      include 'improp.i'
      include 'imptor.i'
      include 'iounit.i'
      include 'korbs.i'
      include 'ktrtor.i'
      include 'kvdws.i'
      include 'math.i'
      include 'mpole3b.i'
      include 'mpole.i'
      include 'mpole3b_3.i'
      include 'mpole3b_2.i'
      include 'opbend.i'
      include 'opdist.i'
      include 'piorbs.i'
      include 'pistuf.i'
      include 'pitors.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'potent.i'
      include 'solute.i'
      include 'strbnd.i'
      include 'strtor.i'
      include 'tors.i'
      include 'tortor.i'
      include 'units.i'
      include 'urey.i'
      include 'vdw.i'
      include 'vdwpot.i'
      integer i,j,k
      integer ia,ib,ic
      integer id,ie,ig
      integer ixaxe
      integer iyaxe
      integer izaxe
      integer fold(6)
      real*8 bla,blc
      real*8 radj,rad4j
      real*8 ampli(6)
      real*8 phase(6)
      real*8 mpl(13)
      logical header
      logical active(*)
c
c
c     number of each type of interaction and site
c
      if (n .ne. 0) then
         write (iout,10)
   10    format (/,' Interactions and Sites :',/)
      end if
      if (nbond .ne. 0) then
         write (iout,20)  nbond
   20    format (' Bond Stretches',19x,i15)
      end if
      if (nangle .ne. 0) then
         write (iout,30)  nangle
   30    format (' Angle Bends',22x,i15)
      end if
      if (nstrbnd .ne. 0) then
         write (iout,40)  nstrbnd
   40    format (' Stretch-Bends',20x,i15)
      end if
      if (nurey .ne. 0) then
         write (iout,50)  nurey
   50    format (' Urey-Bradley',21x,i15)
      end if
      if (nangang .ne. 0) then
         write (iout,50)  nangang
   60    format (' Angle-Angles',21x,i15)
      end if
      if (nopbend .ne. 0) then
         write (iout,70)  nopbend
   70    format (' Out-of-Plane Bends',15x,i15)
      end if
      if (nopdist .ne. 0) then
         write (iout,80)  nopdist
   80    format (' Out-of-Plane Distances',11x,i15)
      end if
      if (niprop .ne. 0) then
         write (iout,90)  niprop
   90    format (' Improper Dihedrals',15x,i15)
      end if
      if (nitors .ne. 0) then
         write (iout,100)  nitors
  100    format (' Improper Torsions',16x,i15)
      end if
      if (ntors .ne. 0) then
         write (iout,110)  ntors
  110    format (' Torsional Angles',17x,i15)
      end if
      if (npitors .ne. 0) then
         write (iout,120)  npitors
  120    format (' Pi-Orbital Torsions',14x,i15)
      end if
      if (nstrtor .ne. 0) then
         write (iout,130)  nstrtor
  130    format (' Stretch-Torsions',17x,i15)
      end if
      if (ntortor .ne. 0) then
         write (iout,140)  ntortor
  140    format (' Torsion-Torsions',17x,i15)
      end if
      if (nvdw .ne. 0) then
         write (iout,150)  nvdw
  150    format (' Van der Waals Sites',14x,i15)
      end if
      if (nion .ne. 0) then
         write (iout,160)  nion
  160    format (' Atomic Partial Charges',11x,i15)
      end if
      if (ndipole .ne. 0) then
         write (iout,170)  ndipole
  170    format (' Bond Dipole Moments',14x,i15)
      end if
      if (npole .ne. 0) then
         write (iout,180)  npole
  180    format (' Atomic Multipoles',16x,i15)
      end if
      if (npolar .ne. 0) then
         write (iout,190)  npolar
  190    format (' Polarizable Sites',16x,i15)
      end if
      if (norbit .ne. 0) then
         write (iout,200)  norbit
  200    format (' Pisystem Atoms',19x,i15)
      end if
      if (nbpi .ne. 0) then
         write (iout,210)  nbpi
  210    format (' Conjugated Pi-Bonds',14x,i15)
      end if
c
c     parameters used for molecular mechanics atom types
c
      header = .true.
      do i = 1, n
         if (active(i)) then
            if (header) then
               header = .false.
               write (iout,220)
  220          format (/,' Atom Type Definition Parameters :',
     &                 //,3x,'Atom',2x,'Symbol',2x,'Type',
     &                    2x,'Class',2x,'Atomic',3x,'Mass',
     &                    2x,'Valence',2x,'Description',/)
            end if
            write (iout,230)  i,name(i),type(i),class(i),atomic(i),
     &                        mass(i),valence(i),story(i)
  230       format (i6,5x,a3,2i7,i6,f10.3,i5,5x,a24)
         end if
      end do
c
c     parameters used for van der Waals interactions
c
      if (use_vdw) then
         header = .true.
         k = 0
         do i = 1, n
            if (active(i)) then
               k = k + 1
               if (header) then
                  header = .false.
                  write (iout,240)
  240             format (/,' Van der Waals Parameters :',
     &                    //,10x,'Atom Number',7x,'Size',
     &                       3x,'Epsilon',3x,'Size 1-4',
     &                       3x,'Eps 1-4',3x,'Reduction',/)
               end if
               j = class(i)
               if (vdwindex .eq. 'TYPE')  j = type(i)
               if (rad(j).eq.rad4(j) .and. eps(j).eq.eps4(j)) then
                  radj = rad(j)
                  if (radsiz .eq. 'DIAMETER')  radj = 2.0d0 * radj
                  if (radtyp .eq. 'SIGMA')  radj = radj / twosix
                  if (reduct(j) .eq. 0.0d0) then
                     write (iout,250)  k,i,radj,eps(j)
  250                format (i6,3x,i6,7x,2f10.4)
                  else
                     write (iout,260)  k,i,radj,eps(j),reduct(j)
  260                format (i6,3x,i6,7x,2f10.4,22x,f10.4)
                  end if
               else
                  radj = rad(j)
                  rad4j = rad4(j)
                  if (radsiz .eq. 'DIAMETER') then
                     radj = 2.0d0 * radj
                     rad4j = 2.0d0 * rad4j
                  end if
                  if (radtyp .eq. 'SIGMA') then
                     radj = radj / twosix
                     rad4j = rad4j / twosix
                  end if
                  if (reduct(j) .eq. 0.0d0) then
                     write (iout,270)  k,i,radj,eps(j),rad4j,eps4(j)
  270                format (i6,3x,i6,7x,2f10.4,1x,2f10.4)
                  else
                     write (iout,280)  k,i,radj,eps(j),rad4j,
     &                                eps4(j),reduct(j)
  280                format (i6,3x,i6,7x,2f10.4,1x,2f10.4,1x,f10.4)
                  end if
               end if
            end if
         end do
      end if
c
c     parameters used for bond stretching interactions
c
      if (use_bond) then
         header = .true.
         do i = 1, nbond
            ia = ibnd(1,i)
            ib = ibnd(2,i)
            if (active(ia) .or. active(ib)) then
               if (header) then
                  header = .false.
                  write (iout,290)
  290             format (/,' Bond Stretching Parameters :',
     &                    //,10x,'Atom Numbers',25x,'KS',7x,'Bond',/)
               end if
               write (iout,300)  i,ia,ib,bk(i),bl(i)
  300          format (i6,3x,2i6,19x,f10.3,f10.4)
            end if
         end do
      end if
c
c     parameters used for angle bending interactions
c
      if (use_angle) then
         header = .true.
         do i = 1, nangle
            ia = iang(1,i)
            ib = iang(2,i)
            ic = iang(3,i)
            if (active(ia) .or. active(ib) .or. active(ic)) then
               if (header) then
                  header = .false.
                  write (iout,310)
  310             format (/,' Angle Bending Parameters :',
     &                    //,13x,'Atom Numbers',22x,'KB',
     &                       6x,'Angle',3x,'Fold',4x,'Type',/)
               end if
               if (angtyp(i) .eq. 'HARMONIC') then
                  write (iout,320)  i,ia,ib,ic,ak(i),anat(i)
  320             format (i6,3x,3i6,13x,2f10.3)
               else if (angtyp(i) .eq. 'IN-PLANE') then
                  write (iout,330)  i,ia,ib,ic,ak(i),anat(i)
  330             format (i6,3x,3i6,13x,2f10.3,9x,'In-Plane')
               else if (angtyp(i) .eq. 'IN-PLANE') then
                  write (iout,340)  i,ia,ib,ic,ak(i),anat(i)
  340             format (i6,3x,3i6,13x,2f10.3,9x,'Linear')
               else if (angtyp(i) .eq. 'FOURIER ') then
                  write (iout,350)  i,ia,ib,ic,ak(i),anat(i),afld(i)
  350             format (i6,3x,3i6,13x,2f10.3,f7.1,2x,'Fourier')
               end if
            end if
         end do
      end if
c
c     parameters used for stretch-bend interactions
c
      if (use_strbnd) then
         header = .true.
         do i = 1, nstrbnd
            k = isb(1,i)
            ia = iang(1,k)
            ib = iang(2,k)
            ic = iang(3,k)
            if (active(ia) .or. active(ib) .or. active(ic)) then
               if (header) then
                  header = .false.
                  write (iout,360)
  360             format (/,' Stretch-Bend Parameters :',
     &                    //,13x,'Atom Numbers',8x,'KSB 1',5x,'KSB 2',
     &                       6x,'Angle',3x,'Bond 1',3x,'Bond 2',/)
               end if
               bla = 0.0d0
               blc = 0.0d0
               if (isb(2,i) .ne. 0)  bla = bl(isb(2,i))
               if (isb(3,i) .ne. 0)  blc = bl(isb(3,i))
               write (iout,370)  i,ia,ib,ic,sbk(1,i),sbk(2,i),
     &                           anat(k),bla,blc
  370          format (i6,3x,3i6,1x,2f10.3,2x,f9.3,2f9.4)
            end if
         end do
      end if
c
c     parameters used for Urey-Bradley interactions
c
      if (use_urey) then
         header = .true.
         do i = 1, nurey
            ia = iury(1,i)
            ib = iury(2,i)
            ic = iury(3,i)
            if (active(ia) .or. active(ic)) then
               if (header) then
                  header = .false.
                  write (iout,380)
  380             format (/,' Urey-Bradley Parameters :',
     &                    //,13x,'Atom Numbers',21x,'KUB',
     &                       4x,'Distance',/)
               end if
               write (iout,390)  i,ia,ib,ic,uk(i),ul(i)
  390          format (i6,3x,3i6,13x,f10.3,f10.4)
            end if
         end do
      end if
c
c     parameters used for out-of-plane bend interactions
c
      if (use_opbend) then
         header = .true.
         do i = 1, nopbend
            k = iopb(i)
            ia = iang(1,k)
            ib = iang(2,k)
            ic = iang(3,k)
            id = iang(4,k)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,400)
  400             format (/,' Out-of-Plane Bend Parameters :',
     &                    //,17x,'Atom Numbers',19x,'KOPB',/)
               end if
               write (iout,410)  i,id,ib,ia,ic,opbk(i)
  410          format (i6,3x,4i6,9x,f10.3)
            end if
         end do
      end if
c
c     parameters used for out-of-plane distance interactions
c
      if (use_opdist) then
         header = .true.
         do i = 1, nopdist
            ia = iopd(1,i)
            ib = iopd(2,i)
            ic = iopd(3,i)
            id = iopd(4,i)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,420)
  420             format (/,' Out-of-Plane Distance Parameters :',
     &                    //,17x,'Atom Numbers',19x,'KOPD',/)
               end if
               write (iout,430)  i,ia,ib,ic,id,opdk(i)
  430          format (i6,3x,4i6,9x,f10.3)
            end if
         end do
      end if
c
c     parameters used for improper dihedral interactions
c
      if (use_improp) then
         header = .true.
         do i = 1, niprop
            ia = iiprop(1,i)
            ib = iiprop(2,i)
            ic = iiprop(3,i)
            id = iiprop(4,i)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,440)
  440             format (/,' Improper Dihedral Parameters :',
     &                    //,17x,'Atom Numbers',19x,'KID',
     &                       4x,'Dihedral',/)
               end if
               write (iout,450)  i,ia,ib,ic,id,kprop(i),vprop(i)
  450          format (i6,3x,4i6,9x,2f10.4)
            end if
         end do
      end if
c
c     parameters used for improper torsion interactions
c
      if (use_imptor) then
         header = .true.
         do i = 1, nitors
            ia = iitors(1,i)
            ib = iitors(2,i)
            ic = iitors(3,i)
            id = iitors(4,i)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,460)
  460             format (/,' Improper Torsion Parameters :',
     &                    //,17x,'Atom Numbers',11x,
     &                       'Amplitude, Phase and Periodicity',/)
               end if
               j = 0
               if (itors1(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 1
                  ampli(j) = itors1(1,i)
                  phase(j) = itors1(2,i)
               end if
               if (itors2(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 2
                  ampli(j) = itors2(1,i)
                  phase(j) = itors2(2,i)
               end if
               if (itors3(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 3
                  ampli(j) = itors3(1,i)
                  phase(j) = itors3(2,i)
               end if
               if (j .eq. 0) then
                  write (iout,470)  i,ia,ib,ic,id
  470             format (i6,3x,4i6)
               else if (j .eq. 1) then
                  write (iout,480)  i,ia,ib,ic,id,
     &                              ampli(1),phase(1),fold(1)
  480             format (i6,3x,4i6,10x,f10.3,f8.1,i4)
               else if (j .eq. 2) then
                  write (iout,490)  i,ia,ib,ic,id,(ampli(k),
     &                              phase(k),fold(k),k=1,j)
  490             format (i6,3x,4i6,2x,2(f10.3,f6.1,i4))
               else
                  write (iout,500)  i,ia,ib,ic,id,(ampli(k),
     &                              nint(phase(k)),fold(k),k=1,j)
  500             format (i6,3x,4i6,4x,3(f8.3,i4,'/',i1))
               end if
            end if
         end do
      end if
c
c     parameters used for torsional interactions
c
      if (use_tors) then
         header = .true.
         do i = 1, ntors
            ia = itors(1,i)
            ib = itors(2,i)
            ic = itors(3,i)
            id = itors(4,i)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,510)
  510             format (/,' Torsional Angle Parameters :',
     &                    //,17x,'Atom Numbers',11x,
     &                       'Amplitude, Phase and Periodicity',/)
               end if
               j = 0
               if (tors1(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 1
                  ampli(j) = tors1(1,i)
                  phase(j) = tors1(2,i)
               end if
               if (tors2(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 2
                  ampli(j) = tors2(1,i)
                  phase(j) = tors2(2,i)
               end if
               if (tors3(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 3
                  ampli(j) = tors3(1,i)
                  phase(j) = tors3(2,i)
               end if
               if (tors4(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 4
                  ampli(j) = tors4(1,i)
                  phase(j) = tors4(2,i)
               end if
               if (tors5(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 5
                  ampli(j) = tors5(1,i)
                  phase(j) = tors5(2,i)
               end if
               if (tors6(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 6
                  ampli(j) = tors6(1,i)
                  phase(j) = tors6(2,i)
               end if
               if (j .eq. 0) then
                  write (iout,520)  i,ia,ib,ic,id
  520             format (i6,3x,4i6)
               else
                  write (iout,530)  i,ia,ib,ic,id,(ampli(k),
     &                              nint(phase(k)),fold(k),k=1,j)
  530             format (i6,3x,4i6,4x,6(f8.3,i4,'/',i1))
               end if
            end if
         end do
      end if
c
c     parameters used for pi-orbital torsion interactions
c
      if (use_pitors) then
         header = .true.
         do i = 1, npitors
            ia = ipit(1,i)
            ib = ipit(2,i)
            ic = ipit(3,i)
            id = ipit(4,i)
            ie = ipit(5,i)
            ig = ipit(6,i)
            if (active(ia) .or. active(ib) .or. active(ic) .or.
     &             active(id) .or. active(ie) .or. active(ig)) then
               if (header) then
                  header = .false.
                  write (iout,540)
  540             format (/,' Pi-Orbital Torsion Parameters :',
     &                    //,10x,'Atom Numbers',19x,'Amplitude',/)
               end if
               write (iout,550)  i,ic,id,kpit(i)
  550          format (i6,3x,2i6,19x,f10.4)
            end if
         end do
      end if
c
c     parameters used for stretch-torsion interactions
c
      if (use_strtor) then
         header = .true.
         do i = 1, nstrtor
            k = ist(1,i)
            ia = itors(1,k)
            ib = itors(2,k)
            ic = itors(3,k)
            id = itors(4,k)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,560)
  560             format (/,' Stretch-Torsion Parameters :',
     &                    //,17x,'Atom Numbers',10x,'Length',
     &                       5x,'Torsion Terms',/)
               end if
               j = 0
               if (kst(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 1
                  ampli(j) = kst(1,i)
                  phase(j) = tors1(2,k)
               end if
               if (kst(2,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 2
                  ampli(j) = kst(2,i)
                  phase(j) = tors2(2,k)
               end if
               if (kst(3,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 3
                  ampli(j) = kst(3,i)
                  phase(j) = tors3(2,k)
               end if
               write (iout,570)  i,ia,ib,ic,id,bl(ist(2,i)),
     &                           (ampli(k),nint(phase(k)),
     &                           fold(k),k=1,j)
  570          format (i6,3x,4i6,2x,f10.4,1x,3(f8.3,i4,'/',i1))
            end if
         end do
      end if
c
c     parameters used for torsion-torsion interactions
c
      if (use_tortor) then
         header = .true.
         do i = 1, ntortor
            k = itt(1,i)
            ia = ibitor(1,k)
            ib = ibitor(2,k)
            ic = ibitor(3,k)
            id = ibitor(4,k)
            ie = ibitor(5,k)
            if (active(ia) .or. active(ib) .or. active(ic) .or.
     &                active(id) .or. active(ie)) then
               if (header) then
                  header = .false.
                  write (iout,580)
  580             format (/,' Torsion-Torsion Parameters :',
     &                    //,20x,'Atom Numbers',18x,'Spline Grid',/)
               end if
               j = itt(2,i)
               write (iout,590)  i,ia,ib,ic,id,ie,tnx(j),tny(j)
  590          format (i6,3x,5i6,10x,2i6)
            end if
         end do
      end if
c
c     parameters used for atomic partial charges
c
      if (use_charge .or. use_chgdpl) then
         header = .true.
         do i = 1, nion
            ia = iion(i)
            ib = jion(i)
            ic = kion(i)
            if (active(ia) .or. active(ic)) then
               if (header) then
                  header = .false.
                  write (iout,600)
  600             format (/,' Atomic Partial Charge Parameters :',
     &                    /,45x,'Neighbor',3x,'Cutoff',
     &                    /,10x,'Atom Number',13x,'Charge',
     &                       7x,'Site',6x,'Site',/)
               end if
               if (ia.eq.ib .and. ia.eq.ic) then
                  write (iout,610)  i,ia,pchg(i)
  610             format (i6,3x,i6,15x,f10.4)
               else
                  write (iout,620)  i,ia,pchg(i),ib,ic
  620             format (i6,3x,i6,15x,f10.4,5x,i6,4x,i6)
               end if
            end if
         end do
      end if
c
c     parameters used for bond dipole moments
c
      if (use_dipole .or. use_chgdpl) then
         header = .true.
         do i = 1, ndipole
            ia = idpl(1,i)
            ib = idpl(2,i)
            if (active(ia) .or. active(ib)) then
               if (header) then
                  header = .false.
                  write (iout,630)
  630             format (/,' Bond Dipole Moment Parameters :',
     &                    //,10x,'Atom Numbers',22x,'Dipole',
     &                       3x,'Position',/)
               end if
               write (iout,640)  i,ia,ib,bdpl(i),sdpl(i)
  640          format (i6,3x,2i6,19x,f10.4,f10.3)
            end if
         end do
      end if
c
c     parameters used for atomic multipole moments
c
      if (use_mpole) then
         header = .true.
         do i = 1, npole
            ia = ipole(i)
            if (active(ia)) then
               if (header) then
                  header = .false.
                  write (iout,650)
  650             format (/,' Atomic Multipole Parameters :',
     &                    //,11x,'Atom',3x,'Z-Axis',1x,'X-Axis',
     &                       1x,'Y-Axis',2x,
     &                       'Frame',11x,'Multipole Moments',/)
               end if
               izaxe = zaxis(i)
               ixaxe = xaxis(i)
               iyaxe = yaxis(i)
               if (iyaxe .lt. 0)  iyaxe = -iyaxe
               mpl(1) = pole(1,i)
               do j = 2, 4
                  mpl(j) = pole(j,i) / bohr
               end do
               do j = 5, 13
                  mpl(j) = 3.0d0 * pole(j,i) / bohr**2
               end do
               if (izaxe .eq. 0) then
                  write (iout,660)  i,ia,polaxe(i),
     &                              (mpl(j),j=1,5),mpl(8),mpl(9),
     &                              (mpl(j),j=11,13)
  660             format (i6,3x,i6,25x,a8,2x,f9.5,/,50x,3f9.5,
     &                    /,50x,f9.5,/,50x,2f9.5,/,50x,3f9.5)
               else if (ixaxe .eq. 0) then
                  write (iout,670)  i,ia,izaxe,polaxe(i),
     &                              (mpl(j),j=1,5),mpl(8),mpl(9),
     &                              (mpl(j),j=11,13)
  670             format (i6,3x,i6,1x,i7,17x,a8,2x,f9.5,/,50x,3f9.5,
     &                    /,50x,f9.5,/,50x,2f9.5,/,50x,3f9.5)
               else  if (iyaxe .eq. 0) then
                  write (iout,680)  i,ia,izaxe,ixaxe,polaxe(i),
     &                              (mpl(j),j=1,5),mpl(8),mpl(9),
     &                              (mpl(j),j=11,13)
  680             format (i6,3x,i6,1x,2i7,10x,a8,2x,f9.5,/,50x,3f9.5,
     &                    /,50x,f9.5,/,50x,2f9.5,/,50x,3f9.5)
               else
                  write (iout,690)  i,ia,izaxe,ixaxe,iyaxe,polaxe(i),
     &                              (mpl(j),j=1,5),mpl(8),mpl(9),
     &                              (mpl(j),j=11,13)
  690             format (i6,3x,i6,1x,3i7,3x,a8,2x,f9.5,/,50x,3f9.5,
     &                    /,50x,f9.5,/,50x,2f9.5,/,50x,3f9.5)
               end if
            end if
         end do
      end if
c
c     parameters used for dipole polarizability
c
      if (use_polar) then
         header = .true.
         do i = 1, npole
            ia = ipole(i)
            if (active(ia)) then
               if (header) then
                  header = .false.
                  write (iout,700)
  700             format (/,' Dipole Polarizability Parameters :',
     &                    //,10x,'Atom Number',5x,'Alpha',5x,'Damp',
     &                       6x,'Polarization Group',/)
               end if
               write (iout,710)  i,ia,polarity(i),thole(i),
     &                           (ip11(j,ia),j=1,np11(ia))
  710          format (i6,3x,i6,6x,f10.4,f9.3,3x,20i6)
            end if
         end do
      end if
c
c     parameters used for empirical solvation
c
      if (use_solv) then
         header = .true.
         k = 0
         do i = 1, n
            if (active(i)) then
               k = k + 1
               if (header) then
                  header = .false.
                  write (iout,720)
  720             format (/,' Empirical Solvation Parameters :',
     &                    //,10x,'Atom Number',13x,'Radius',
     &                       3x,'ASP Value',/)
               end if
               write (iout,730)  k,i,rsolv(i),asolv(i)
  730          format (i6,3x,i6,15x,2f10.4)
            end if
         end do
      end if
c
c     parameters used for conjugated pisystem atoms
c
      if (use_orbit) then
         header = .true.
         do i = 1, norbit
            ia = iorbit(i)
            j = class(ia)
            if (header) then
               header = .false.
               write (iout,740)
  740          format (/,' Conjugated Pi-Atom Parameters :',
     &                 //,10x,'Atom Number',14x,'Nelect',
     &                    6x,'Ionize',4x,'Repulsion',/)
            end if
            write (iout,750)  i,ia,electron(j),ionize(j),repulse(j)
  750       format (i6,3x,i6,17x,f8.1,3x,f10.4,2x,f10.4)
         end do
      end if
c
c     parameters used for conjugated pibond interactions
c
      if (use_orbit) then
         header = .true.
         do i = 1, nbpi
            ia = ibpi(2,i)
            ib = ibpi(3,i)
            if (header) then
               header = .false.
               write (iout,760)
  760          format (/,' Conjugated Pi-Bond Parameters :',
     &                 //,10x,'Atom Numbers',21x,'K Slope',
     &                    3x,'L Slope',/)
            end if
            write (iout,770)  i,ia,ib,kslope(i),lslope(i)
  770       format (i6,3x,2i6,19x,2f10.4)
         end do
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine enrgyze  --  compute & report energy analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "enrgyze" is an auxiliary routine for the analyze program
c     that performs the energy analysis and prints the total and
c     intermolecular energies
c
c
      subroutine enrgyze ()
      implicit none
      include 'sizes3b.i'
      include 'sizes.i'
      include 'atoms.i'
      include 'atoms3b.i'
      include 'cutoff.i'
      include 'inform.i'
      include 'inter.i'
      include 'iounit.i'
      include 'molcul.i'
      include 'molcul3b.i'
      include 'combo.i'
      integer i,j,epvr1,epvr2
      real*8 xr,yr,zr,r2
      real*8 energy
      character*120 fstr

c
c     perform the energy analysis by atom and component
c
      call analysis (energy)
c      call findr2
c      moli = 1
c      do moli1 = 1, nmol 
c         do moli2 = moli1+1, nmol 
c            if (r2b(moli) .ge. 5.0d0) then
c               print*,r2b(moli)
c               write(epvr1,*)r2b(moli),ep2analyze(moli)
c                write(epvr1,*)ep2analyze(moli)
c            end if
c            moli = moli + 1
c         end do
c      end do
c      call findr3
c      moli = 1
c      do moli1 = 1, nmol 
c         do moli2 = moli1+1, nmol
c            do moli3 = moli2+1, nmol
c               r2=r3b_1(moli)+r3b_2(moli)+r3b_3(moli) 
c               write(epvr2,*)r2,ep3analyze(moli)
c                write(epvr2,*)ep3analyze(moli)
c               moli = moli + 1
c            end do
c         end do
c       end do
c
c
c     print out the total potential energy of the system
c
      fstr = '(/,'' Total Potential Energy :'',8x,f16.4,'' Kcal/mole'')'
      if (digits .ge. 6)  fstr(32:39) = '6x,f18.6'
      if (digits .ge. 8)  fstr(32:39) = '4x,f20.8'
      if (abs(energy) .ge. 1.0d10)  fstr(35:35) = 'd'
      write (iout,fstr)  energy
c
c     intermolecular energy for systems with multiple molecules
c
      fstr = '(/,'' Intermolecular Energy :'',9x,f16.4,'' Kcal/mole'')'
      if (digits .ge. 6)  fstr(31:38) = '7x,f18.6'
      if (digits .ge. 8)  fstr(31:38) = '5x,f20.8'
      if (abs(einter) .ge. 1.0d10)  fstr(34:34) = 'd'
      if (nmol.gt.1 .and. nmol.lt.n .and. .not.use_ewald)
     &   write (iout,fstr)  einter
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine partyze  --  energy component decomposition  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "partyze" prints the energy component and number of
c     interactions for each of the potential energy terms
c
c
      subroutine partyze
      include 'action.i'
      include 'cutoff.i'
      include 'energi.i'
      include 'inform.i'
      include 'iounit.i'
      include 'potent.i'
      character*12 form1
      character*12 form2
      character*120 fstr
c
c
c     write out each energy component to the desired precision
c
      form1 = '5x,f16.4,i15'
      if (digits .ge. 6)  form1 = '3x,f18.6,i15'
      if (digits .ge. 8)  form1 = '1x,f20.8,i15'
      form2 = form1(1:3)//'d'//form1(5:12)
      fstr = '(/,'' Energy Component Breakdown :'',
     &          11x,''Kcal/mole'',6x,''Interactions''/)'
      write (iout,fstr)
      if (use_bond .and. neb.ne.0) then
         fstr = '('' Bond Stretching'',12x,'//form1//')'
         write (iout,fstr)  eb,neb
      end if
      if (use_angle .and. nea.ne.0) then
         fstr = '('' Angle Bending'',14x,'//form1//')'
         write (iout,fstr)  ea,nea
      end if
      if (use_strbnd .and. neba.ne.0) then
         fstr = '('' Stretch-Bend'',15x,'//form1//')'
         write (iout,fstr)  eba,neba
      end if
      if (use_urey .and. neub.ne.0) then
         fstr = '('' Urey-Bradley'',15x,'//form1//')'
         write (iout,fstr)  eub,neub
      end if
      if (use_angang .and. neaa.ne.0) then
         fstr = '('' Angle-Angle'',16x,'//form1//')'
         write (iout,fstr)  eaa,neaa
      end if
      if (use_opbend .and. neopb.ne.0) then
         fstr = '('' Out-of-Plane Bend'',10x,'//form1//')'
         write (iout,fstr)  eopb,neopb
      end if
      if (use_opdist .and. neopd.ne.0) then
         fstr = '('' Out-of-Plane Distance'',6x,'//form1//')'
         write (iout,fstr)  eopd,neopd
      end if
      if (use_improp .and. neid.ne.0) then
         fstr = '('' Improper Dihedral'',10x,'//form1//')'
         write (iout,fstr)  eid,neid
      end if
      if (use_imptor .and. neit.ne.0) then
         fstr = '('' Improper Torsion'',11x,'//form1//')'
         write (iout,fstr)  eit,neit
      end if
      if (use_tors .and. net.ne.0) then
         fstr = '('' Torsional Angle'',12x,'//form1//')'
         write (iout,fstr)  et,net
      end if
      if (use_pitors .and. nept.ne.0) then
         fstr = '('' Pi-Orbital Torsion'',9x,'//form1//')'
         write (iout,fstr)  ept,nept
      end if
      if (use_strtor .and. nebt.ne.0) then
         fstr = '('' Stretch-Torsion'',12x,'//form1//')'
         write (iout,fstr)  ebt,nebt
      end if
      if (use_tortor .and. nett.ne.0) then
         fstr = '('' Torsion-Torsion'',12x,'//form1//')'
         write (iout,fstr)  ett,nett
      end if
      if (use_vdw .and. nev.ne.0) then
         if (abs(ev) .lt. 1.0d10) then
            fstr = '('' Van der Waals'',14x,'//form1//')'
         else
            fstr = '('' Van der Waals'',14x,'//form2//')'
         end if
         write (iout,fstr)  ev,nev
      end if
      if (use_charge .and. nec.ne.0) then
         if (abs(ec) .lt. 1.0d10) then
            fstr = '('' Charge-Charge'',14x,'//form1//')'
         else
            fstr = '('' Charge-Charge'',14x,'//form2//')'
         end if
         write (iout,fstr)  ec,nec
      end if
      if (use_chgdpl .and. necd.ne.0) then
         if (abs(ecd) .lt. 1.0d10) then
            fstr = '('' Charge-Dipole'',14x,'//form1//')'
         else
            fstr = '('' Charge-Dipole'',14x,'//form2//')'
         end if
         write (iout,fstr)  ecd,necd
      end if
      if (use_dipole .and. ned.ne.0) then
         if (abs(ed) .lt. 1.0d10) then
            fstr = '('' Dipole-Dipole'',14x,'//form1//')'
         else
            fstr = '('' Dipole-Dipole'',14x,'//form2//')'
         end if
         write (iout,fstr)  ed,ned
      end if
      if (use_mpole .and. nem.ne.0) then
         if (abs(em) .lt. 1.0d10) then
            fstr = '('' Atomic Multipoles'',10x,'//form1//')'
         else
            fstr = '('' Atomic Multipoles'',10x,'//form2//')'
         end if
         write (iout,fstr)  em,nem
      end if
      if (use_polar .and. nep.ne.0) then
         if (abs(ep) .lt. 1.0d10) then
            fstr = '('' Polarization'',15x,'//form1//')'
         else
            fstr = '('' Polarization'',15x,'//form2//')'
         end if
         write (iout,fstr)  ep,nep
      end if
      if (use_rxnfld .and. ner.ne.0) then
         fstr = '('' Reaction Field'',13x,'//form1//')'
         write (iout,fstr)  er,ner
      end if
      if (use_solv .and. nes.ne.0) then
         fstr = '('' Implicit Solvation'',9x,'//form1//')'
         write (iout,fstr)  es,nes
      end if
      if (use_metal .and. nelf.ne.0) then
         fstr = '('' Metal Ligand Field'',9x,'//form1//')'
         write (iout,fstr)  elf,nelf
      end if
      if (use_geom .and. neg.ne.0) then
         fstr = '('' Geometric Restraints'',7x,'//form1//')'
         write (iout,fstr)  eg,neg
      end if
      if (use_extra .and. nex.ne.0) then
         fstr = '('' Extra Energy Terms'',9x,'//form1//')'
         write (iout,fstr)  ex,nex
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine propyze  --  electrostatic & inertial analysis  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "propyze" finds and prints the total charge, dipole moment
c     components, radius of gyration, moments of inertia, internal
c     virial and pressure
c
c
      subroutine propyze
      include 'sizes3b.i'
      include 'sizes.i'
      include 'atoms.i'
      include 'atoms3b.i'
      include 'chgpot.i'
      include 'chgpot3b.i'
      include 'iounit.i'
      include 'moment.i'
      include 'virial.i'
      real*8 rg,energy
      real*8, allocatable :: derivs(:,:)
c
c
c     get the total charge, dipole and quadrupole moments
c
      call moments
      write (iout,10)  netchg
   10 format (/,' Total Electric Charge :',13x,f12.5,' Electrons')
      write (iout,20)  netdpl,xdpl,ydpl,zdpl
   20 format (/,' Dipole Moment Magnitude :',11x,f12.3,' Debyes',
     &        //,' Dipole X,Y,Z-Components :',11x,3f12.3)
      write (iout,30)  xxqdp,xyqdp,xzqdp,yxqdp,yyqdp,
     &                 yzqdp,zxqdp,zyqdp,zzqdp
   30 format (/,' Quadrupole Moment Tensor :',10x,3f12.3,
     &        /,6x,'(Buckinghams)',18x,3f12.3,
     &        /,37x,3f12.3)
      write (iout,40)  netqdp(1),netqdp(2),netqdp(3)
   40 format (/,' Principal Axes Quadrupole :',9x,3f12.3)
      if (dielec .ne. 1.0d0) then
         write (iout,50)  dielec
   50    format (/,' Dielectric Constant :',15x,f12.3)
         write (iout,60)  netchg/sqrt(dielec)
   60    format (' Effective Total Charge :',12x,f12.5,' Electrons')
         write (iout,70)  netdpl/sqrt(dielec)
   70    format (' Effective Dipole Moment :',11x,f12.3,' Debyes')
      end if
c
c     get the radius of gyration and moments of inertia
c
      call gyrate (rg)
      write (iout,80)  rg
   80 format (/,' Radius of Gyration :',16x,f12.3,' Angstroms')
      call inertia (1)
c
c     get the internal virial tensor via gradient calculation
c
      allocate (derivs(3,n))
      call gradient (energy,derivs)
      deallocate (derivs)
      write (iout,90)  (vir(1,i),vir(2,i),vir(3,i),i=1,3)
   90 format (/,' Internal Virial Tensor :',12x,3f12.3,
     &        /,37x,3f12.3,/,37x,3f12.3)
c
c     get two alternative dE/dV values and a pressure estimate
c
      call ptest
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine atomyze  --  individual atom energy analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "atomyze" prints the potential energy components broken
c     down by atom and to a choice of precision
c
c
      subroutine atomyze (active)
      include 'sizes3b.i'
      include 'sizes.i'
      include 'analyz.i'
      include 'atoms.i'
      include 'atoms3b.i'
      include 'inform.i'
      include 'iounit.i'
      integer i
      logical active(*)
      character*120 fstr
c
c
c     energy partitioning over the individual atoms
c
      fstr = '(/,'' Potential Energy Breakdown over Atoms :'')'
      write (iout,fstr)
      if (digits .ge. 8) then
         write (iout,10)
   10    format (/,'  Atom',9x,'EB',14x,'EA',14x,'EBA',13x,'EUB',
     &           /,15x,'EAA',13x,'EOPB',12x,'EOPD',12x,'EID',
     &           /,15x,'EIT',13x,'ET',14x,'EPT',13x,'EBT',
     &           /,15x,'ETT',13x,'EV',14x,'EC',14x,'ECD',
     &           /,15x,'ED',14x,'EM',14x,'EP',14x,'ER',
     &           /,15x,'ES',14x,'ELF',13x,'EG',14x,'EX')
      else if (digits .ge. 6) then
         write (iout,20)
   20    format (/,'  Atom',8x,'EB',12x,'EA',12x,'EBA',11x,'EUB',
     &              11x,'EAA',
     &           /,14x,'EOPB',10x,'EOPD',10x,'EID',11x,'EIT',
     &              11x,'ET',
     &           /,14x,'EPT',11x,'EBT',11x,'ETT',11x,'EV',12x,'EC',
     &           /,14x,'ECD',11x,'ED',12x,'EM',12x,'EP',12x,'ER',
     &           /,14x,'ES',12x,'ELF',11x,'EG',12x,'EX')
      else
         write (iout,30)
   30    format (/,'  Atom',8x,'EB',10x,'EA',10x,'EBA',9x,'EUB',
     &              9x,'EAA',9x,'EOPB',
     &           /,14x,'EOPD',8x,'EID',9x,'EIT',9x,'ET',10x,'EPT',
     &              9x,'EBT',
     &           /,14x,'ETT',9x,'EV',10x,'EC',10x,'ECD',9x,'ED',
     &              10x,'EM',
     &           /,14x,'EP',10x,'ER',10x,'ES',10x,'ELF',9x,'EG',
     &              10x,'EX')
      end if
      if (digits .ge. 8) then
         fstr = '(/,i6,4f16.8,/,6x,4f16.8,/,6x,4f16.8,'//
     &             '/,6x,4f16.8,/,6x,4f16.8,/,6x,4f16.8)'
      else if (digits .ge. 6) then
         fstr = '(/,i6,5f14.6,/,6x,5f14.6,/,6x,5f14.6,'//
     &             '/,6x,5f14.6,/,6x,4f14.6)'
      else
         fstr = '(/,i6,6f12.4,/,6x,6f12.4,/,6x,6f12.4,/,6x,6f12.4)'
      end if
      do i = 1, n
         if (active(i)) then
            write (iout,fstr)  i,aeb(i),aea(i),aeba(i),aeub(i),aeaa(i),
     &                         aeopb(i),aeopd(i),aeid(i),aeit(i),aet(i),
     &                         aept(i),aebt(i),aett(i),aev(i),aec(i),
     &                         aecd(i),aed(i),aem(i),aep(i),aer(i),
     &                         aes(i),aelf(i),aeg(i),aex(i)
         end if
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine amberyze  --  parameter format for Amber setup  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "amberyze" prints the force field parameters in a format needed
c     by the Amber setup protocol for using AMOEBA within Amber; this
c     is essentially the "paramyze" format from TINKER 4.3
c
c
      subroutine amberyze (active)
      implicit none
      include 'sizes3b.i'
      include 'sizes.i'
      include 'angang.i'
      include 'angle.i'
      include 'angpot.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'atoms3b.i'
      include 'bitor.i'
      include 'bond.i'
      include 'charge.i'
      include 'dipole.i'
      include 'improp.i'
      include 'imptor.i'
      include 'iounit.i'
      include 'korbs.i'
      include 'ktrtor.i'
      include 'kvdws.i'
      include 'math.i'      
      include 'mpole3b.i'
      include 'mpole.i'
      include 'mpole3b_3.i'
      include 'mpole3b_2.i'
      include 'opbend.i'
      include 'opdist.i'
      include 'piorbs.i'
      include 'pistuf.i'
      include 'pitors.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'potent.i'
      include 'solute.i'
      include 'strbnd.i'
      include 'strtor.i'
      include 'tors.i'
      include 'tortor.i'
      include 'units.i'
      include 'urey.i'
      include 'vdw.i'
      include 'vdwpot.i'
      integer i,j,k
      integer ia,ib,ic
      integer id,ie,ig
      integer ixaxe
      integer iyaxe
      integer izaxe
      integer fold(6)
      real*8 bla,blc
      real*8 sbavg
      real*8 radj,rad4j
      real*8 ampli(6)
      real*8 phase(6)
      real*8 mpl(13)
      logical header
      logical active(*)
c
c
c     number of each type of interaction and site
c
      if (n .ne. 0) then
         write (iout,10)
   10    format (/,' Total Numbers of Atoms and Interactions :')
         write (iout,20)  n
   20    format (/,' Atoms in System',11x,i15)
      end if
      if (norbit .ne. 0) then
         write (iout,30)  norbit
   30    format (' Pisystem Atoms',12x,i15)
      end if
      if (nbond .ne. 0) then
         write (iout,40)  nbond
   40    format (' Bond Stretches',12x,i15)
      end if
      if (nbpi .ne. 0) then
         write (iout,50)  nbpi
   50    format (' Conjugated Pi-Bonds',7x,i15)
      end if
      if (nangle .ne. 0) then
         write (iout,60)  nangle
   60    format (' Angle Bends',15x,i15)
      end if
      if (nstrbnd .ne. 0) then
         write (iout,70)  nstrbnd
   70    format (' Stretch-Bends',13x,i15)
      end if
      if (nurey .ne. 0) then
         write (iout,80)  nurey
   80    format (' Urey-Bradley',14x,i15)
      end if
      if (nangang .ne. 0) then
         write (iout,80)  nangang
   90    format (' Angle-Angles',14x,i15)
      end if
      if (nopbend .ne. 0) then
         write (iout,100)  nopbend
  100    format (' Out-of-Plane Bends',8x,i15)
      end if
      if (nopdist .ne. 0) then
         write (iout,110)  nopdist
  110    format (' Out-of-Plane Distances',4x,i15)
      end if
      if (niprop .ne. 0) then
         write (iout,120)  niprop
  120    format (' Improper Dihedrals',8x,i15)
      end if
      if (nitors .ne. 0) then
         write (iout,130)  nitors
  130    format (' Improper Torsions',9x,i15)
      end if
      if (ntors .ne. 0) then
         write (iout,140)  ntors
  140    format (' Torsional Angles',10x,i15)
      end if
      if (npitors .ne. 0) then
         write (iout,150)  npitors
  150    format (' Pi-Orbital Torsions',7x,i15)
      end if
      if (nstrtor .ne. 0) then
         write (iout,160)  nstrtor
  160    format (' Stretch-Torsions',10x,i15)
      end if
      if (ntortor .ne. 0) then
         write (iout,170)  ntortor
  170    format (' Torsion-Torsions',10x,i15)
      end if
      if (nion .ne. 0) then
         write (iout,180)  nion
  180    format (' Atomic Partial Charges',4x,i15)
      end if
      if (ndipole .ne. 0) then
         write (iout,190)  ndipole
  190    format (' Bond Dipole Moments',7x,i15)
      end if
      if (npole .ne. 0) then
         write (iout,200)  npole
  200    format (' Polarizable Multipoles',4x,i15)
      end if
c
c     parameters used for molecular mechanics atom types
c
      header = .true.
      do i = 1, n
         if (active(i)) then
            if (header) then
               header = .false.
               write (iout,220)
  220          format (/,' Atom Type Definition Parameters :',
     &                 //,3x,'Atom',2x,'Symbol',2x,'Type',
     &                    2x,'Class',2x,'Atomic',3x,'Mass',
     &                    2x,'Valence',2x,'Description',/)
            end if
            write (iout,230)  i,name(i),type(i),class(i),atomic(i),
     &                        mass(i),valence(i),story(i)
  230       format (i6,5x,a3,2i7,i6,f10.3,i5,5x,a24)
         end if
      end do
c
c     parameters used for van der Waals interactions
c
      if (use_vdw) then
         header = .true.
         k = 0
         do i = 1, n
            if (active(i)) then
               k = k + 1
               if (header) then
                  header = .false.
                  write (iout,240)
  240             format (/,' Van der Waals Parameters :',
     &                    //,10x,'Atom Number',7x,'Radius',
     &                       3x,'Epsilon',3x,'Rad 1-4',
     &                       3x,'Eps 1-4',3x,'Reduction',/)
               end if
               j = class(i)
               if (vdwindex .eq. 'TYPE')  j = type(i)
               if (rad(j).eq.rad4(j) .and. eps(j).eq.eps4(j)) then
                  radj = rad(j)
                  if (radtyp .eq. 'SIGMA')  radj = radj / twosix
                  if (reduct(j) .eq. 0.0d0) then
                     write (iout,250)  k,i,radj,eps(j)
  250                format (i6,3x,i6,7x,2f10.4)
                  else
                     write (iout,260)  k,i,radj,eps(j),reduct(j)
  260                format (i6,3x,i6,7x,2f10.4,22x,f10.4)
                  end if
               else
                  radj = rad(j)
                  rad4j = rad4(j)
                  if (radtyp .eq. 'SIGMA') then
                     radj = radj / twosix
                     rad4j = rad4j / twosix
                  end if
                  if (reduct(j) .eq. 0.0d0) then
                     write (iout,270)  k,i,radj,eps(j),rad4j,eps4(j)
  270                format (i6,3x,i6,7x,2f10.4,1x,2f10.4)
                  else
                     write (iout,280)  k,i,radj,eps(j),rad4j,
     &                                eps4(j),reduct(j)
  280                format (i6,3x,i6,7x,2f10.4,1x,2f10.4,1x,f10.4)
                  end if
               end if
            end if
         end do
      end if
c
c     parameters used for bond stretching interactions
c
      if (use_bond) then
         header = .true.
         do i = 1, nbond
            ia = ibnd(1,i)
            ib = ibnd(2,i)
            if (active(ia) .or. active(ib)) then
               if (header) then
                  header = .false.
                  write (iout,290)
  290             format (/,' Bond Stretching Parameters :',
     &                    //,10x,'Atom Numbers',25x,'KS',7x,'Length',/)
               end if
               write (iout,300)  i,ia,ib,bk(i),bl(i)
  300          format (i6,3x,2i6,19x,f10.3,f10.4)
            end if
         end do
      end if
c
c     parameters used for angle bending interactions
c
      if (use_angle) then
         header = .true.
         do i = 1, nangle
            ia = iang(1,i)
            ib = iang(2,i)
            ic = iang(3,i)
            if (active(ia) .or. active(ib) .or. active(ic)) then
               if (header) then
                  header = .false.
                  write (iout,310)
  310             format (/,' Angle Bending Parameters :',
     &                    //,13x,'Atom Numbers',22x,'KB',
     &                       6x,'Angle',3x,'Fold',4x,'Type',/)
               end if
               if (angtyp(i) .eq. 'HARMONIC') then
                  write (iout,320)  i,ia,ib,ic,ak(i),anat(i)
  320             format (i6,3x,3i6,13x,2f10.3)
               else if (angtyp(i) .eq. 'IN-PLANE') then
                  write (iout,330)  i,ia,ib,ic,ak(i),anat(i)
  330             format (i6,3x,3i6,13x,2f10.3,9x,'In-Plane')
               else if (angtyp(i) .eq. 'IN-PLANE') then
                  write (iout,340)  i,ia,ib,ic,ak(i),anat(i)
  340             format (i6,3x,3i6,13x,2f10.3,9x,'Linear')
               else if (angtyp(i) .eq. 'FOURIER ') then
                  write (iout,350)  i,ia,ib,ic,ak(i),anat(i),afld(i)
  350             format (i6,3x,3i6,13x,2f10.3,f7.1,2x,'Fourier')
               end if
            end if
         end do
      end if
c
c     parameters used for stretch-bend interactions
c
      if (use_strbnd) then
         header = .true.
         do i = 1, nstrbnd
            k = isb(1,i)
            ia = iang(1,k)
            ib = iang(2,k)
            ic = iang(3,k)
            if (active(ia) .or. active(ib) .or. active(ic)) then
               if (header) then
                  header = .false.
                  write (iout,360)
  360             format (/,' Stretch-Bend Parameters :',
     &                    //,13x,'Atom Numbers',11x,'KSB',
     &                       6x,'Angle',3x,'Length1',3x,'Length2',/)
               end if
               bla = 0.0d0
               blc = 0.0d0
               if (isb(2,i) .ne. 0)  bla = bl(isb(2,i))
               if (isb(3,i) .ne. 0)  blc = bl(isb(3,i))
               sbavg = (sbk(1,i)+sbk(2,i)) * 0.5d0
               write (iout,370)  i,ia,ib,ic,sbavg,anat(k),bla,blc
  370          format (i6,3x,3i6,f13.4,3f10.4)
            end if
         end do
      end if
c
c     parameters used for Urey-Bradley interactions
c
      if (use_urey) then
         header = .true.
         do i = 1, nurey
            ia = iury(1,i)
            ib = iury(2,i)
            ic = iury(3,i)
            if (active(ia) .or. active(ic)) then
               if (header) then
                  header = .false.
                  write (iout,380)
  380             format (/,' Urey-Bradley Parameters :',
     &                    //,13x,'Atom Numbers',21x,'KUB',
     &                       4x,'Distance',/)
               end if
               write (iout,390)  i,ia,ib,ic,uk(i),ul(i)
  390          format (i6,3x,3i6,13x,f10.3,f10.4)
            end if
         end do
      end if
c
c     parameters used for out-of-plane bend interactions
c
      if (use_opbend) then
         header = .true.
         do i = 1, nopbend
            k = iopb(i)
            ia = iang(1,k)
            ib = iang(2,k)
            ic = iang(3,k)
            id = iang(4,k)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,400)
  400             format (/,' Out-of-Plane Bending Parameters :',
     &                    //,17x,'Atom Numbers',19x,'KOPB',/)
               end if
               opbk(i) = opbk(i) * (opbunit/.02191418)
               write (iout,410)  i,id,ib,ia,ic,opbk(i)
  410          format (i6,3x,4i6,9x,f10.4)
            end if
         end do
      end if
c
c     parameters used for out-of-plane distance interactions
c
      if (use_opdist) then
         header = .true.
         do i = 1, nopdist
            ia = iopd(1,i)
            ib = iopd(2,i)
            ic = iopd(3,i)
            id = iopd(4,i)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,420)
  420             format (/,' Out-of-Plane Distance Parameters :',
     &                    //,17x,'Atom Numbers',19x,'KOPD',/)
               end if
               write (iout,430)  i,ia,ib,ic,id,opdk(i)
  430          format (i6,3x,4i6,9x,f10.3)
            end if
         end do
      end if
c
c     parameters used for improper dihedral interactions
c
      if (use_improp) then
         header = .true.
         do i = 1, niprop
            ia = iiprop(1,i)
            ib = iiprop(2,i)
            ic = iiprop(3,i)
            id = iiprop(4,i)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,440)
  440             format (/,' Improper Dihedral Parameters :',
     &                    //,17x,'Atom Numbers',19x,'KID',
     &                       4x,'Dihedral',/)
               end if
               write (iout,450)  i,ia,ib,ic,id,kprop(i),vprop(i)
  450          format (i6,3x,4i6,9x,2f10.4)
            end if
         end do
      end if
c
c     parameters used for improper torsion interactions
c
      if (use_imptor) then
         header = .true.
         do i = 1, nitors
            ia = iitors(1,i)
            ib = iitors(2,i)
            ic = iitors(3,i)
            id = iitors(4,i)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,460)
  460             format (/,' Improper Torsion Parameters :',
     &                    //,17x,'Atom Numbers',11x,
     &                       'Amplitude, Phase and Periodicity',/)
               end if
               j = 0
               if (itors1(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 1
                  ampli(j) = itors1(1,i)
                  phase(j) = itors1(2,i)
               end if
               if (itors2(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 2
                  ampli(j) = itors2(1,i)
                  phase(j) = itors2(2,i)
               end if
               if (itors3(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 3
                  ampli(j) = itors3(1,i)
                  phase(j) = itors3(2,i)
               end if
               if (j .eq. 0) then
                  write (iout,470)  i,ia,ib,ic,id
  470             format (i6,3x,4i6)
               else if (j .eq. 1) then
                  write (iout,480)  i,ia,ib,ic,id,
     &                              ampli(1),phase(1),fold(1)
  480             format (i6,3x,4i6,10x,f10.3,f8.1,i4)
               else if (j .eq. 2) then
                  write (iout,490)  i,ia,ib,ic,id,(ampli(k),
     &                              phase(k),fold(k),k=1,j)
  490             format (i6,3x,4i6,2x,2(f10.3,f6.1,i4))
               else
                  write (iout,500)  i,ia,ib,ic,id,(ampli(k),
     &                              nint(phase(k)),fold(k),k=1,j)
  500             format (i6,3x,4i6,4x,3(f8.3,i4,'/',i1))
               end if
            end if
         end do
      end if
c
c     parameters used for torsional interactions
c
      if (use_tors) then
         header = .true.
         do i = 1, ntors
            ia = itors(1,i)
            ib = itors(2,i)
            ic = itors(3,i)
            id = itors(4,i)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,510)
  510             format (/,' Torsional Angle Parameters :',
     &                    //,17x,'Atom Numbers',11x,
     &                       'Amplitude, Phase and Periodicity',/)
               end if
               j = 0
               if (tors1(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 1
                  ampli(j) = tors1(1,i)
                  phase(j) = tors1(2,i)
               end if
               if (tors2(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 2
                  ampli(j) = tors2(1,i)
                  phase(j) = tors2(2,i)
               end if
               if (tors3(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 3
                  ampli(j) = tors3(1,i)
                  phase(j) = tors3(2,i)
               end if
               if (tors4(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 4
                  ampli(j) = tors4(1,i)
                  phase(j) = tors4(2,i)
               end if
               if (tors5(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 5
                  ampli(j) = tors5(1,i)
                  phase(j) = tors5(2,i)
               end if
               if (tors6(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 6
                  ampli(j) = tors6(1,i)
                  phase(j) = tors6(2,i)
               end if
               if (j .eq. 0) then
                  write (iout,520)  i,ia,ib,ic,id
  520             format (i6,3x,4i6)
               else
                  write (iout,530)  i,ia,ib,ic,id,(ampli(k),
     &                              nint(phase(k)),fold(k),k=1,j)
  530             format (i6,3x,4i6,4x,6(f8.3,i4,'/',i1))
               end if
            end if
         end do
      end if
c
c     parameters used for pi-orbital torsion interactions
c
      if (use_pitors) then
         header = .true.
         do i = 1, npitors
            ia = ipit(1,i)
            ib = ipit(2,i)
            ic = ipit(3,i)
            id = ipit(4,i)
            ie = ipit(5,i)
            ig = ipit(6,i)
            if (active(ia) .or. active(ib) .or. active(ic) .or.
     &             active(id) .or. active(ie) .or. active(ig)) then
               if (header) then
                  header = .false.
                  write (iout,540)
  540             format (/,' Pi-Orbital Torsion Parameters :',
     &                    //,10x,'Atom Numbers',19x,'Amplitude',/)
               end if
               write (iout,550)  i,ic,id,kpit(i)
  550          format (i6,3x,2i6,19x,f10.4)
            end if
         end do
      end if
c
c     parameters used for stretch-torsion interactions
c
      if (use_strtor) then
         header = .true.
         do i = 1, nstrtor
            k = ist(1,i)
            ia = itors(1,k)
            ib = itors(2,k)
            ic = itors(3,k)
            id = itors(4,k)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,560)
  560             format (/,' Stretch-Torsion Parameters :',
     &                    //,17x,'Atom Numbers',10x,'Length',
     &                       5x,'Torsion Terms',/)
               end if
               j = 0
               if (kst(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 1
                  ampli(j) = kst(1,i)
                  phase(j) = tors1(2,k)
               end if
               if (kst(2,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 2
                  ampli(j) = kst(2,i)
                  phase(j) = tors2(2,k)
               end if
               if (kst(3,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 3
                  ampli(j) = kst(3,i)
                  phase(j) = tors3(2,k)
               end if
               write (iout,570)  i,ia,ib,ic,id,bl(ist(2,i)),
     &                           (ampli(k),nint(phase(k)),
     &                           fold(k),k=1,j)
  570          format (i6,3x,4i6,2x,f10.4,1x,3(f8.3,i4,'/',i1))
            end if
         end do
      end if
c
c     parameters used for torsion-torsion interactions
c
      if (use_tortor) then
         header = .true.
         do i = 1, ntortor
            k = itt(1,i)
            ia = ibitor(1,k)
            ib = ibitor(2,k)
            ic = ibitor(3,k)
            id = ibitor(4,k)
            ie = ibitor(5,k)
            if (active(ia) .or. active(ib) .or. active(ic) .or.
     &                active(id) .or. active(ie)) then
               if (header) then
                  header = .false.
                  write (iout,580)
  580             format (/,' Torsion-Torsion Parameters :',
     &                    //,20x,'Atom Numbers',18x,'Spline Grid',/)
               end if
               j = itt(2,i)
               write (iout,590)  i,ia,ib,ic,id,ie,tnx(j),tny(j)
  590          format (i6,3x,5i6,10x,2i6)
            end if
         end do
      end if
c
c     parameters used for atomic partial charges
c
      if (use_charge .or. use_chgdpl) then
         header = .true.
         do i = 1, nion
            ia = iion(i)
            ib = jion(i)
            ic = kion(i)
            if (active(ia) .or. active(ic)) then
               if (header) then
                  header = .false.
                  write (iout,600)
  600             format (/,' Atomic Partial Charge Parameters :',
     &                    /,45x,'Neighbor',3x,'Cutoff',
     &                    /,10x,'Atom Number',13x,'Charge',
     &                       7x,'Site',6x,'Site',/)
               end if
               if (ia.eq.ib .and. ia.eq.ic) then
                  write (iout,610)  i,ia,pchg(i)
  610             format (i6,3x,i6,15x,f10.4)
               else
                  write (iout,620)  i,ia,pchg(i),ib,ic
  620             format (i6,3x,i6,15x,f10.4,5x,i6,4x,i6)
               end if
            end if
         end do
      end if
c
c     parameters used for bond dipole moments
c
      if (use_dipole .or. use_chgdpl) then
         header = .true.
         do i = 1, ndipole
            ia = idpl(1,i)
            ib = idpl(2,i)
            if (active(ia) .or. active(ib)) then
               if (header) then
                  header = .false.
                  write (iout,630)
  630             format (/,' Bond Dipole Moment Parameters :',
     &                    //,10x,'Atom Numbers',22x,'Dipole',
     &                       3x,'Position',/)
               end if
               write (iout,640)  i,ia,ib,bdpl(i),sdpl(i)
  640          format (i6,3x,2i6,19x,f10.4,f10.3)
            end if
         end do
      end if
c
c     parameters used for atomic multipole moments
c
      if (use_mpole) then
         header = .true.
         do i = 1, npole
            ia = ipole(i)
            if (active(ia)) then
               if (header) then
                  header = .false.
                  write (iout,650)
  650             format (/,' Atomic Multipole Parameters :',
     &                    //,12x,'Atom',4x,'Coordinate Frame',
     &                       ' Definition',7x,'Multipole Moments',/)
               end if
               izaxe = zaxis(i)
               ixaxe = xaxis(i)
               iyaxe = yaxis(i)
               if (iyaxe .lt. 0)  iyaxe = -iyaxe
               mpl(1) = pole(1,i)
               do j = 2, 4
                  mpl(j) = pole(j,i) / bohr
               end do
               do j = 5, 13
                  mpl(j) = 3.0d0 * pole(j,i) / bohr**2
               end do
               if (izaxe .eq. 0) then
                  write (iout,660)  i,ia,polaxe(i),
     &                              (mpl(j),j=1,5),mpl(8),mpl(9),
     &                              (mpl(j),j=11,13)
  660             format (i6,3x,i6,25x,a8,2x,f9.5,/,50x,3f9.5,
     &                    /,50x,f9.5,/,50x,2f9.5,/,50x,3f9.5)
               else if (ixaxe .eq. 0) then
                  write (iout,670)  i,ia,izaxe,polaxe(i),
     &                              (mpl(j),j=1,5),mpl(8),mpl(9),
     &                              (mpl(j),j=11,13)
  670             format (i6,3x,i6,1x,i7,17x,a8,2x,f9.5,/,50x,3f9.5,
     &                    /,50x,f9.5,/,50x,2f9.5,/,50x,3f9.5)
               else  if (iyaxe .eq. 0) then
                  write (iout,680)  i,ia,izaxe,ixaxe,polaxe(i),
     &                              (mpl(j),j=1,5),mpl(8),mpl(9),
     &                              (mpl(j),j=11,13)
  680             format (i6,3x,i6,1x,2i7,10x,a8,2x,f9.5,/,50x,3f9.5,
     &                    /,50x,f9.5,/,50x,2f9.5,/,50x,3f9.5)
               else
                  write (iout,690)  i,ia,izaxe,ixaxe,iyaxe,polaxe(i),
     &                              (mpl(j),j=1,5),mpl(8),mpl(9),
     &                              (mpl(j),j=11,13)
  690             format (i6,3x,i6,1x,3i7,3x,a8,2x,f9.5,/,50x,3f9.5,
     &                    /,50x,f9.5,/,50x,2f9.5,/,50x,3f9.5)
               end if
            end if
         end do
      end if
c
c     parameters used for dipole polarizability
c
c
      if (use_polar) then
         header = .true.
         do i = 1, npole
            ia = ipole(i)
            if (active(ia)) then
               if (header) then
                  header = .false.
                  write (iout,700)
  700             format (/,' Dipole Polarizability Parameters :',
     &                    //,10x,'Atom Number',9x,'Alpha',8x,
     &                       'Polarization Group',/)
               end if
               write (iout,710)  i,ia,polarity(i),
     &                           (ip11(j,ia),j=1,np11(ia))
  710          format (i6,3x,i6,10x,f10.4,5x,20i6)
            end if
         end do
      end if
c
c     parameters used for empirical solvation
c
      if (use_solv) then
         header = .true.
         k = 0
         do i = 1, n
            if (active(i)) then
               k = k + 1
               if (header) then
                  header = .false.
                  write (iout,720)
  720             format (/,' Empirical Solvation Parameters :',
     &                    //,10x,'Atom Number',13x,'Radius',
     &                       3x,'ASP Value',/)
               end if
               write (iout,730)  k,i,rsolv(i),asolv(i)
  730          format (i6,3x,i6,15x,2f10.4)
            end if
         end do
      end if
c
c     parameters used for conjugated pisystem atoms
c
      if (use_orbit) then
         header = .true.
         do i = 1, norbit
            ia = iorbit(i)
            j = class(ia)
            if (header) then
               header = .false.
               write (iout,740)
  740          format (/,' Conjugated Pi-Atom Parameters :',
     &                 //,10x,'Atom Number',14x,'Nelect',
     &                    6x,'Ionize',4x,'Repulsion',/)
            end if
            write (iout,750)  i,ia,electron(j),ionize(j),repulse(j)
  750       format (i6,3x,i6,17x,f8.1,3x,f10.4,2x,f10.4)
         end do
      end if
c
c     parameters used for conjugated pibond interactions
c
      if (use_orbit) then
         header = .true.
         do i = 1, nbpi
            ia = ibpi(2,i)
            ib = ibpi(3,i)
            if (header) then
               header = .false.
               write (iout,760)
  760          format (/,' Conjugated Pi-Bond Parameters :',
     &                 //,10x,'Atom Numbers',21x,'K Slope',
     &                    3x,'L Slope',/)
            end if
            write (iout,770)  i,ia,ib,kslope(i),lslope(i)
  770       format (i6,3x,2i6,19x,2f10.4)
         end do
      end if
      return
      end
