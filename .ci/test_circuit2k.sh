export BASE=$DATADIR/cellular/circuit-2k

srun appositionizer \
    --output=$PWD \
    --appositionspace=axodendritic \
    --recipe $BASE/bioname/recipe.json \
    --circuit-config $BASE/circuit_config.json

set +x
hash module 2> /dev/null || . /etc/profile.d/modules.sh
. /gpfs/bbp.cscs.ch/apps/hpc/jenkins/config/modules.sh
module load unstable py-pytouchreader
module list
set -xeo pipefail
touchdiff --full --segment-precision 1 --distance-margin 1.0 --distance-precision=0.1 {$BASE/touches/binary,.}/touches.0 | head -n 500
