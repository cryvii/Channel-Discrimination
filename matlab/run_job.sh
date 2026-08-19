#!/bin/bash
#SBATCH --job-name=upper_bound_test
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8000
#SBATCH --time=06:00:00
#SBATCH --output=slurm_logs/%x_%j.out
#SBATCH --error=slurm_logs/%x_%j.err

set -uo pipefail

MATLAB_BIN="/home/u224567/MATLAB/bin/matlab"
SCRIPT_DIR="/home/u224567/MATLAB/SDP/3channel-final/matlab/"   # adjust if the pipeline lives elsewhere
REPO_ROOT="/home/u224567/MATLAB/SDP/3channel-final"   # adjust if the pipeline lives elsewhere
RESULTS_DIR="${REPO_ROOT}/results"
RAW_DIR="${RESULTS_DIR}/raw"
PROC_DIR="${RESULTS_DIR}/processed"
PYTHON_BIN="python3"

# ---------------------------------------------------------------------------
# RUN_ID is computed HERE (not inside MATLAB) so both the MATLAB save path
# and the Python post-processing step agree on the same folder name, keyed
# by the SLURM job ID. This is what makes "job ID + N + sample" a valid
# lookup key for find_channel.py.
# ---------------------------------------------------------------------------
TAG="${SLURM_JOB_NAME:-job}"
RUN_ID="$(date +%Y%m%d_%H%M%S)_job${SLURM_JOB_ID:-local}_${TAG}"

mkdir -p "${RAW_DIR}" "${PROC_DIR}" "$(dirname "${0}")/slurm_logs" 2>/dev/null

echo "============================================================"
echo "  Job: ${SLURM_JOB_NAME}  (ID: ${SLURM_JOB_ID})"
echo "  Node: $(hostname)"
echo "  Started: $(date)"
echo "  N: ${N:-<unset>}"
echo "  RUN_ID: ${RUN_ID}"
echo "  Raw dir: ${RAW_DIR}/${RUN_ID}"
echo "============================================================"

"${MATLAB_BIN}" -nodisplay -nosplash -nodesktop -r "addpath(genpath('/home/u224567/MATLAB/SDP')); addpath(genpath('${SCRIPT_DIR}')); cd('${SCRIPT_DIR}'); generate_results('protocols',[1 2 5 6 8 9],'N_start',${N},'N_max',${N},'n_samples',50,'seed','random','channels','diff','channel_type','channels','savedir','${RAW_DIR}','run_id','${RUN_ID}'); exit"

MATLAB_EXIT=$?
echo "MATLAB finished: $(date)  |  exit code: ${MATLAB_EXIT}"

# Move the SLURM logs into the run folder now that we know it exists,
# instead of leaving them scattered in the repo root.
if [ -d "${RAW_DIR}/${RUN_ID}" ]; then
    [ -f "slurm_logs/${TAG}_${SLURM_JOB_ID}.out" ] && cp "slurm_logs/${TAG}_${SLURM_JOB_ID}.out" "${RAW_DIR}/${RUN_ID}/slurm.out"
    [ -f "slurm_logs/${TAG}_${SLURM_JOB_ID}.err" ] && cp "slurm_logs/${TAG}_${SLURM_JOB_ID}.err" "${RAW_DIR}/${RUN_ID}/slurm.err"
fi

if [ ${MATLAB_EXIT} -ne 0 ]; then
    echo "MATLAB failed (exit ${MATLAB_EXIT}) -- skipping Python post-processing."
    exit ${MATLAB_EXIT}
fi

# ---------------------------------------------------------------------------
# Automatic post-processing: raw .mat -> tidy CSVs + plots + dashboard.
# Zero manual steps for a submitted job to produce clean, reproducible results.
# ---------------------------------------------------------------------------
if [ -d "${RAW_DIR}/${RUN_ID}" ]; then
    echo "Running Python processor on ${RAW_DIR}/${RUN_ID} ..."
    "${PYTHON_BIN}" "${REPO_ROOT}/python/process_run.py" \
        --run-dir "${RAW_DIR}/${RUN_ID}" \
        --output-dir "${PROC_DIR}/${RUN_ID}"
    PROC_EXIT=$?
    echo "Python processor finished: $(date)  |  exit code: ${PROC_EXIT}"
else
    echo "WARNING: expected raw run dir ${RAW_DIR}/${RUN_ID} not found; skipping processing."
    PROC_EXIT=1
fi

echo "Finished: $(date)"
exit $(( MATLAB_EXIT != 0 ? MATLAB_EXIT : PROC_EXIT ))