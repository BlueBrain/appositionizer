#!/bin/bash
#SBATCH --time=01:30:00
#SBATCH --account=proj16
#SBATCH --constraint=cpu
#SBATCH --partition=prod
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=40

export APPO_CIRCUIT="S1.v6a"
export APPO_MODE="modern"

source "${APPO_DIR}/sbatch/run_appositionizer.sh"
