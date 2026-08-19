#!/usr/bin/env bash
#
# run_pipeline.sh
#
# Runs the whole data pipeline end-to-end, in dependency order:
#
#   1. merge_samples.py   matlab/slurm_logs/*.out  -> Data/unitary/, Data/channels/
#   2. process_data.py    Data/unitary/, Data/channels/ -> results/processed/<class>/
#                            long_format.csv, pivot_mean_p.csv, pivot_mean_t.csv,
#                            pivot_sample_count.csv, pivot_solved_rate.csv
#   3. analyze_data.py    results/processed/<class>/ -> ratio.csv, max_separation.csv
#
# Usage (run from anywhere -- it locates the repo root itself):
#   ./run_pipeline.sh
#
# This script lives in python/, next to merge_samples.py, process_data.py,
# and analyze_data.py.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=============================================="
echo " 3channel-final data pipeline"
echo " Working directory: $SCRIPT_DIR"
echo "=============================================="

echo
echo "----- [1/3] merge_samples.py -----"
python3 merge_samples.py

echo
echo "----- [2/3] process_data.py -----"
python3 process_data.py

echo
echo "----- [3/3] analyze_data.py -----"
python3 analyze_data.py

echo
echo "=============================================="
echo " Done. Outputs written to:"
echo "   ../Data/unitary/                (per-N merged .txt, from merge_samples.py)"
echo "   ../Data/channels/               (per-N merged .txt, from merge_samples.py)"
echo "   ../results/processed/unitary/   (long_format.csv, pivot_*.csv, ratio.csv, max_separation.csv)"
echo "   ../results/processed/channels/  (same, for the 'channels' data)"
echo "=============================================="