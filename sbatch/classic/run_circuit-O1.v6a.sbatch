#!/bin/bash
#SBATCH --time=01:30:00
#SBATCH --account=proj16
#SBATCH --constraint=cpu
#SBATCH --partition=prod
#SBATCH --ntasks=320

export APPO_CIRCUIT="O1.v6a"
export APPO_MODE="classic"

source "${APPO_DIR}/sbatch/run_appositionizer.sh"
