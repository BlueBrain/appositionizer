export BASE=$DATADIR/cellular/circuit-1k

echo "SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MPI_OPENMP_INTEROP=Yes

srun appositionizer \
    --tmp-dir=$PWD/tmp_dir \
    --output=$PWD \
    --appositionspace=axodendritic \
    --recipe $BASE/bioname/recipe.json \
    --circuit-config $BASE/circuit_config_merged.json

set +x
hash module 2> /dev/null || . /etc/profile.d/modules.sh
. /gpfs/bbp.cscs.ch/apps/hpc/jenkins/config/modules.sh
module load unstable py-pytouchreader
module list
set -xeo pipefail
touchdiff --full --segment-precision 1 --distance-margin 1.0 --distance-precision=0.1 {$BASE/touches/binary,.}/touches.0 | head -n 500
