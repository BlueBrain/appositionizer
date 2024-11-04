#!/bin/bash
#SBATCH --account=proj16
#SBATCH --constraint=cpu
#SBATCH --partition=prod
#SBATCH --ntasks=32

export APPO_CIRCUIT="1k"
export APPO_MODE="classic"

source "${APPO_DIR}/sbatch/run_appositionizer.sh"
