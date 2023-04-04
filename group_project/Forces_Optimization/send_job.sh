#!/usr/bin/env bash


#SBATCH --job-name=F_opt_ljmd       		# The name of your job, you'll se it in squeue.
#SBATCH --mail-type=NONE              	# Mail events (NONE, BEGIN, END, FAIL, ALL). Sends you an email when the job begins, ends, or fails; you can combine options.
#SBATCH --mail-user=mgisonni@sissa.it   # Where to send the mail

#SBATCH --ntasks=1                   # Number of MPI ranks (1 for MPI serial job)
#SBATCH --cpus-per-task=40            # Number of threads per MPI rank (MAX: 2x32 cores on _partition_2, 2x20 cores on _partition_1) 
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks-per-node=1          # How many tasks on each node
#SBATCH --ntasks-per-socket=1        # How many tasks on each socket
#SBATCH --ntasks-per-core=1          # How many tasks on each core (set to 1 to be sure that different tasks run on different cores on multi-threaded systems)

#SBATCH --mem=900mb                 # Memory per node (MAX: 63500 on the new ones, 40000 on the old ones); incompatible with --mem-per-cpu.


#SBATCH --partition=regular1         # Partition (queue). Avail: regular1, regular2, long1, long2, wide1, wide2, gpu1, gpu2. Multiple partitions are possible.
#SBATCH --time=00:05:00              # Time limit hrs:min:sec
#SBATCH --output=%x.o%j              # Standard output log in TORQUE-style -- WARNING: %x requires a new enough SLURM. Use %j for regular jobs and %A-%a for array jobs
#SBATCH --error=%x.e%j               # Standard error  log in TORQUE-style -- WARNING: %x requires a new enough SLURM. Use %j for regular jobs and %A-%a for array jobs



# ==== Modules part (load all the modules) ===== #

	module load cmake/3.15.4
	module load intel/18.0.3.222

# ==== End of Modules part (load all the modules) ===== #



# ==== Info part (say things) ===== #
#
NOW=`date +%H:%M-%a-%d/%b/%Y`
echo '------------------------------------------------------'
echo 'This job is allocated on '$SLURM_JOB_CPUS_PER_NODE' cpu(s)'
echo 'Job is running on node(s): '
echo  $SLURM_JOB_NODELIST
echo '------------------------------------------------------'
echo 'WORKINFO:'
echo 'SLURM: job starting at           '$NOW
echo 'SLURM: sbatch is running on      '$SLURM_SUBMIT_HOST
echo 'SLURM: executing on cluster      '$SLURM_CLUSTER_NAME
echo 'SLURM: executing on partition    '$SLURM_JOB_PARTITION
echo 'SLURM: working directory is      '$SLURM_SUBMIT_DIR
echo 'SLURM: current home directory is '$(getent passwd $SLURM_JOB_ACCOUNT | cut -d: -f6)
echo ""
echo 'JOBINFO:'
echo 'SLURM: job identifier is         '$SLURM_JOBID
echo 'SLURM: job name is               '$SLURM_JOB_NAME
echo ""
echo 'NODEINFO:'
echo 'SLURM: number of nodes is        '$SLURM_JOB_NUM_NODES
echo 'SLURM: number of cpus/node is    '$SLURM_JOB_CPUS_PER_NODE
echo 'SLURM: number of gpus/node is    '$SLURM_GPUS_PER_NODE
echo '------------------------------------------------------'
#
# ==== End of Info part (say things) ===== #


# # Should not be necessary anymore with SLURM, as this is the default, but you never know...
# cd $SLURM_SUBMIT_DIR


# ==== JOB COMMANDS ===== #

rm -rf build
mkdir build
cd build
cmake ..
cmake --build .
cd examples
make VERBOSE=1;
./../ljmd-serial.x < argon_108.inp 
./../ljmd-serial.x < argon_2916.inp 

# ==== END OF JOB COMMANDS ===== #



# Wait for processes, if any.
echo "Waiting for all the processes to finish..."
wait
