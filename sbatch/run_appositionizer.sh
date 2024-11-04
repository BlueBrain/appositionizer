if [[ -z ${APPO_DIR} ]]
then
    echo "APPO_DIR not set."
    exit -1
fi

# Will load additional env variables, e.g. for Caliper.
if [[ $# -eq 1 ]]
then
    source "$1"
fi


source "${APPO_DIR}/sbatch/${APPO_MODE}/env.sh"
source "${APPO_DIR}/sbatch/circuit-${APPO_CIRCUIT}.sh"

APPO_OUTPUT_DIR=${APPO_OUTPUT_DIR:-${PWD}}
APPO_JOBSIZE=${APPO_JOBSIZE:-4}

time srun appositionizer \
    --output=${APPO_OUTPUT_DIR} \
    --jobsize=${APPO_JOBSIZE} \
    ${APPO_MODERN} \
    --from ${APPO_FROM} \
    ${APPO_FROM_NODESET_FLAG} \
    --to ${APPO_TO} \
    ${APPO_TO_NODESET_FLAG} \
    --appositionspace=${APPO_SPACE} \
    ${APPO_RECIPE} \
    ${APPO_MORPHOLOGIES}
