
# The "typical" command line arguments for Appositionizer are listed in this file. If a
# particular circuit needs something else it can first source this file and then
# overwrite variables as needed.
#
# Typical means:
#  - axodendritic
#  - stored in "the" CI folder
#  - reads morphologies from 'morphologies/h5'
#  - output to a scratch directory

export APPO_BASE=/gpfs/bbp.cscs.ch/project/proj12/connectome-establishment/${APPO_CIRCUIT}
export APPO_FROM="${APPO_BASE}/nodes.h5"
export APPO_TO="${APPO_BASE}/nodes.h5"
export APPO_SPACE=axodendritic 
export APPO_RECIPE="${APPO_BASE}/bioname/builderRecipeAllPathways.xml"
export APPO_MORPHOLOGIES="${APPO_BASE}/morphologies/h5"

if [[ -z ${APPO_OUTPUT_DIR} ]]
then
    SCRATCH="/gpfs/bbp.cscs.ch/project/proj12/scratch/connectome-establishment"
    mkdir -p "${SCRATCH}/appo_output"
    export APPO_OUTPUT_DIR=$(mktemp -d ${SCRATCH}/appo_output/${APPO_CIRCUIT}-XXXXXX)
fi
