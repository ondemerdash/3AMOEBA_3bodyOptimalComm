#PBS -q debug
#PBS -l nodes=2:ppn=8
#PBS -N small16
#PBS -j oe
#PBS -o small16.$PBS_JOBID.log

module swap pgi intel
module swap openmpi openmpi-intel

cd $PBS_O_WORKDIR
mpirun ../analyze_parallel_bcast_evenload_bsend_pme_parall.x -k 512_cluster_no_openmp.key 512_cluster.xyz E > small16.$PBS_JOBID.out
