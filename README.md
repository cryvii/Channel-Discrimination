# Semidefinite Programming for Multi-Slot Quantum Channel Discrimination

This repository contains the numerical infrastructure for evaluating quantum channel discrimination under different causal structures [1]. It provides a complete pipeline to formulate semidefinite programs (SDPs) for nine distinct causal classes, solve them at scale on a SLURM-managed cluster, and process, analyze, and visualize the resulting datasets [1].

The implementation is designed for a fixed number of slots $k=3$ and qubit systems $d=2$, meaning each process matrix and tester lives in $\mathcal{L}(\mathcal{H}^{\otimes 6})$ with dimension $d^{2k}=64$ [1].

---

## Repository Structure

```text
├── channel_bench.html            # Standalone browser-based plotting tool (JS + Chart.js)
├── matlab/
│   ├── channel_discrimination_3copies_primal_9classes.m  # Core SDP dispatch function
│   ├── generate_results.m        # Main orchestration script for channel generation & solving
│   ├── draw_one_channel.m        # Random channel generator (unitary, general, or Pauli)
│   ├── submit_all.sh             # SLURM array job submission script
│   ├── run_job.sh                # SLURM batch execution script
│   ├── TR.m                      # Helper for trace-and-replace operations
│   └── slurm_logs/               # Directory for raw cluster execution logs
├── python/
│   ├── run_pipeline.sh           # Main orchestration script for post-processing
│   ├── merge_samples.py          # Log parser and sample deduplicator (Stage 3)
│   ├── process_data.py           # Data aggregator and pivot table generator (Stage 4)
│   └── analyze_data.py           # Statistical analysis (ratios and maximum separations)
├── Data/                         # Deduplicated text datasets grouped by N
│   ├── unitary/
│   └── channels/
└── results/
    ├── raw/                      # Raw simulation parameters, Choi matrices, and run metadata
    └── processed/                # Aggregated CSVs, pivot tables, and analysis outputs
```

---

## Theoretical Background: Causal Classes

The implementation supports nine protocols/causal classes from the channel discrimination hierarchy [1]:

| Index | Class | Description |
| :---: | :--- | :--- |
| **1** | `PAR` | Parallel strategy |
| **2** | `SEQ` | Sequential strategy (comb-constrained SDP solved over $3!=6$ permutations) |
| **3** | `QC-convFO` | Convex mixtures of fixed-order sequential strategies |
| **4** | `QC-NICC` | Quantum circuits with non-interactive classical control |
| **5** | `pQC-CC` | Purified quantum circuits with classical control |
| **6** | `QC-SupFO` | Quantum circuits with superposition of input-output orders |
| **7** | `QC-NIQC` | Quantum circuits with non-interactive quantum control |
| **8** | `QC-QC` | General quantum circuits with quantum control |
| **9** | `GEN` | General process matrices (unconstrained except for valid testers) |

---

## Pipeline Overview

The data generation and analysis workflow is split into five modular stages [1].

```
[Stage 1: Configure] ──> [Stage 2: Compute] ──> [Stage 3: Merge] ──> [Stage 4: Process] ──> [Stage 5: Visualize]
```

### Stage 1: Configure
The configuration for execution is managed via the SLURM batch script `matlab/run_job.sh` [1]. Within this file, you can customize the parameters passed to `generate_results.m` [1]:
*   Target protocols (e.g., `[1 2 5 6 8 9]`)
*   Number of samples per $N$ (`n_samples`)
*   Channel type (`'unitary'`, `'channels'`/`'superoperator'`, or `'pauli'`)
*   Discrimination mode (`'diff'` for independent Haar-drawn channels, or `'same'`/`'identical'` for a single re-used channel)
*   Random seed

### Stage 2: Compute (Cluster Execution)
Run the automated submission script specifying the lower and upper bounds for the number of hypotheses ($N$) [1]:
```bash
cd matlab
./submit_all.sh 2 35
```
This spawns one independent SLURM array job per value of $N$ [1]. 
*   Each job requests 4 CPUs, 8 GB of RAM, and a 6-hour walltime [1].
*   Outputs are saved to a run-specific directory `results/raw/<RUN_ID>/` containing the exact Choi matrices used (`channels_used_diff.mat/.txt`), raw solver results (`results_diff.txt`), and execution metadata (`meta.json`) [1].
*   Solver progress is outputted to `matlab/slurm_logs/<jobname>_<jobid>.out` [1].

*Note: The script `run_job.sh` attempts to call `python/process_run.py` automatically upon completion [1]. Because this script is omitted from the distribution, the automatic per-job trigger will fail silently [1]. Post-processing should instead be run collectively via Stage 4 [1].*

### Stage 3: Merge and Deduplicate
To parse the raw output files and extract valid samples, run the Python merge utility [1]:
```bash
python python/merge_samples.py
```
This script scans all log files in `matlab/slurm_logs/`, matches completed sample blocks, groups them by $N$ and channel type, and writes consolidated datasets to `Data/unitary/<N>_<n>.txt` or `Data/channels/<N>_<n>.txt` [1].
*   Each block is annotated with a provenance line containing its original `source`, `job_id`, `run_id`, and original sample index to prevent data loss or duplicate counting [1].
*   Stale files from previous runs are cleaned up automatically [1].

### Stage 4: Process and Analyze
Run the central post-processing bash wrapper [1]:
```bash
./python/run_pipeline.sh
```
This script runs the following steps sequentially [1]:
1.  Executes `merge_samples.py` (Stage 3).
2.  Runs `process_data.py` to create a master database (`long_format.csv`) containing every individual sample and its provenance, alongside pivoted summary tables (`pivot_mean_p.csv`, `pivot_mean_t.csv`, `pivot_sample_count.csv`, `pivot_solved_rate.csv`).
3.  Runs `analyze_data.py` to compute pairwise success probability ratios (`ratio.csv`) and find the individual samples producing the largest absolute separation between any two protocols (`max_separation.csv`).

All resulting files are written directly to `results/processed/unitary/` or `results/processed/channels/` [1].

### Stage 5: Visualize
The repository includes a zero-dependency, browser-based plotting tool, `channel_bench.html` (built with vanilla JS and Chart.js) [1].
1.  Open `channel_bench.html` in any web browser [1].
2.  Drag and drop any of the processed pivot CSV files (e.g., `pivot_mean_p.csv`) or `long_format.csv` directly into the window [1].
3.  Use the interface to choose the independent variable, toggle specific protocol traces, select logarithmic scaling, apply line-smoothing, or adjust aesthetics for colorblind-friendly academic publishing [1].
4.  Export vector graphics directly as `.svg` files [1].

---

## Data Traceability

Every final result in `results/processed/` is fully auditable [1]. The pipeline preserves data lineage at every step [1]:
$$\text{Pivot Table / Analysis File} \longrightarrow \text{long\_format.csv} \longrightarrow \text{Data/<class>/<N>\_<n>.txt} \longrightarrow \text{slurm.out} \longrightarrow \text{results/raw/<RUN\_ID>/}$$

By looking up the provenance comment of any sample inside `Data/`, you can retrieve the exact random seed, physical parameters, and generated Choi matrices (`channels_used_diff.mat`) that produced the data point [1].

---

## Dependencies

*   **MATLAB** (R2020a or later recommended)
*   **CVX** (configured with the **Mosek** solver [1])
*   **QETLAB** (for `PartialTrace` and `PermuteSystems` functions [1])
*   **Python 3.x** (with `numpy` and `pandas` libraries [1])
