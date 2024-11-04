#!/bin/bash
#SBATCH --account=proj16
#SBATCH --constraint=cpu
#SBATCH --partition=prod
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=16

export APPO_CIRCUIT="2k"
export APPO_MODE="modern"

source "${APPO_DIR}/sbatch/run_appositionizer.sh"
