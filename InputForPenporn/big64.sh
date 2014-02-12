#PBS -q debug
#PBS -l nodes=8:ppn=8
#PBS -N big64
#PBS -j oe
#PBS -o big64.$PBS_JOBID.log

module swap pgi intel
module swap openmpi openmpi-intel

cd $PBS_O_WORKDIR
mpirun ../analyze_parallel_bcast_evenload_bsend_pme_parall.x -k 6400_cluster_no_openmp.key 6400_cluster.xyz E > big64.$PBS_JOBID.out
