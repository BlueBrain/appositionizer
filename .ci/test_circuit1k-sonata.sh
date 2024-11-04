export BASE=$DATADIR/cellular/circuit-1k

srun appositionizer \
    --output=$PWD \
    --output-population=All \
    --appositionspace=axodendritic \
    --recipe $BASE/bioname/recipe.json \
    --circuit-config $BASE/circuit_config.json

set +x
hash module 2> /dev/null || . /etc/profile.d/modules.sh
. /gpfs/bbp.cscs.ch/apps/hpc/jenkins/config/modules.sh
module load unstable py-functionalizer parquetk
module list
set -xeo pipefail
python ${SPACK_SOURCE_DIR}/tools/sonata_compare.py $BASE/touches/td_edges.h5 edges.h5
