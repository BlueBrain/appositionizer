#!/bin/bash
#SBATCH --time=01:30:00
#SBATCH --account=proj16
#SBATCH --constraint=cpu
#SBATCH --partition=prod
#SBATCH --ntasks=320

# This circuit likely needs longer to run.
# Please update accordingly.

export APPO_CIRCUIT="S1.v6a"
export APPO_MODE="classic"

source "${APPO_DIR}/sbatch/run_appositionizer.sh"
