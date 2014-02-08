#!/bin/bash

#$ -S /bin/bash
#$ -N MPI6400
#$ -q all.q
#$ -pe orte 64 
#$ -cwd

module load ompi-intel 
mpirun /home/ondemerdash/GoodVersionsNoEwaldThgLoop/DontUseReplica_UseImagesNoCutOrigLoopNoFuncNoEwald3bNo1bodyOptEwaldPermOMPParallelTiming/MPIVersionNewest3List_CobarNonDiscrete8/OptSmallList/FewerArrayAlloc/ParallelizeAll/analyze_parallel_bcast_evenload_bsend_pme_parall.x -k /home/ondemerdash/Xyz_and_Keyfiles/6400_cluster_no_openmp.key /home/ondemerdash/Xyz_and_Keyfiles/6400_cluster.xyz E > ParAllPMEFewerArrayAllocBsendOptSmallListlog6400_DontUseReplica_UseImagesNoCutOrigLoopNoFuncNoEwald3bNo1bodyOptEwaldPermOMPParallelTiming_MPIVersionNewest3List_CobarNonDiscrete8_64thread 
