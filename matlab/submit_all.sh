#!/bin/bash
# submit_all.sh — submits run_job.sh once per N value, N=2..35 by default.
# Usage: ./submit_all.sh [N_LOW] [N_HIGH]

N_LOW="${1:-2}"
N_HIGH="${2:-35}"

mkdir -p slurm_logs

for N in $(seq "${N_LOW}" "${N_HIGH}"); do
    sbatch --job-name="unitary_upper_bound_N${N}" \
           --output="slurm_logs/unitary_upper_bound_N${N}_%j.out" \
           --error="slurm_logs/unitary_upper_bound_N${N}_%j.err" \
           --export=ALL,N=${N} \
           run_job.sh
done
