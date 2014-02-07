c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  group3b.i  --  partitioning of system into atom groups  ##
c     ##                                                        ##
c     ############################################################
c
c
c     use_group   flag to use partitioning of system into groups
c     use_intra   flag to include only intragroup interactions
c     use_inter   flag to include only intergroup interactions
c
c
      logical use_group
      logical use_intra
      logical use_inter
      common /group3b/ use_group,use_intra,use_inter
