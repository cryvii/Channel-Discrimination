# Quantum Channel Discrimination under Causal Constraints

This repository contains the numerical infrastructure for evaluating quantum channel discrimination under different causal structures. It provides a complete pipeline to formulate semidefinite programs (SDPs) for nine distinct causal classes, solve them at scale on a SLURM-managed cluster, and process, analyze, and visualize the resulting datasets.

The implementation is designed for a fixed number of slots $k=3$ and qubit systems $d=2$, meaning each process matrix and tester lives in $\mathcal{L}(\mathcal{H}^{\otimes 6})$ with dimension $d^{2k}=64$.

---

## Repository Structure

    ├── channel\_bench.html            # Standalone browser-based plotting tool (JS + Chart.js)
    ├── matlab/
    │   ├── channel\_discrimination\_3copies\_primal\_9classes.m  # Core SDP dispatch function
    │   ├── generate\_results.m        # Main orchestration script for channel generation & solving
    │   ├── draw\_one\_channel.m        # Random channel generator (unitary, general, or Pauli)
    │   ├── submit\_all.sh             # SLURM array job submission script
    │   ├── run\_job.sh                # SLURM batch execution script
    │   ├── TR.m                      # Helper for trace-and-replace operations
    │   └── slurm\_logs/               # Directory for raw cluster execution logs
    ├── python/
    │   ├── run\_pipeline.sh           # Main orchestration script for post-processing
    │   ├── merge\_samples.py          # Log parser and sample deduplicator (Stage 3)
    │   ├── process\_data.py           # Data aggregator and pivot table generator (Stage 4)
    │   └── analyze\_data.py           # Statistical analysis (ratios and maximum separations)
    ├── Data/                         # Deduplicated text datasets grouped by N
    │   ├── unitary/
    │   └── channels/
    └── results/
        ├── raw/                      # Raw simulation parameters, Choi matrices, and run metadata
        └── processed/                # Aggregated CSVs, pivot tables, and analysis outputs

---

## Theoretical Background: Causal Classes

The implementation supports nine protocols/causal classes from the channel discrimination hierarchy:

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

The data generation and analysis workflow is split into five modular stages.

[Stage 1: Configure] ──> [Stage 2: Compute] ──> [Stage 3: Merge] ──> [Stage 4: Process] ──> [Stage 5: Visualize]

### Stage 1: Configure
The configuration for execution is managed via the SLURM batch script `matlab/run\_job.sh`. Within this file, you can customize the parameters passed to `generate\_results.m` [1]:
*   Target protocols (e.g., `[1 2 5 6 8 9]`)
*   Number of samples per $N$ (`n\_samples`)
*   Channel type (`'unitary'`, `'channels'`/`'superoperator'`, or `'pauli'`)
*   Discrimination mode (`'diff'` for independent Haar-drawn channels, or `'same'`/`'identical'` for a single re-used channel)
*   Random seed

### Stage 2: Compute (Cluster Execution)
Run the automated submission script specifying the lower and upper bounds for the number of hypotheses ($N$):

    cd matlab
    ./submit\_all.sh 2 35

This spawns one independent SLURM array job per value of $N$. 
*   Each job requests 4 CPUs, 8 GB of RAM, and a 6-hour walltime.
*   Outputs are saved to a run-specific directory `results/raw/<RUN\_ID>/` containing the exact Choi matrices used (`channels\_used\_diff.mat/.txt`), raw solver results (`results\_diff.txt`), and execution metadata (`meta.json`).
*   Solver progress is outputted to `matlab/slurm\_logs/<jobname>\_<jobid>.out`.

*Note: The script `run\_job.sh` attempts to call `python/process\_run.py` automatically upon completion [1]. Because this script is omitted from the distribution, the automatic per-job trigger will fail silently. Post-processing should instead be run collectively via Stage 4.*

### Stage 3: Merge and Deduplicate
To parse the raw output files and extract valid samples, run the Python merge utility:

    python python/merge\_samples.py

This script scans all log files in `matlab/slurm\_logs/`, matches completed sample blocks, groups them by $N$ and channel type, and writes consolidated datasets to `Data/unitary/<N>\_<n>.txt` or `Data/channels/<N>\_<n>.txt`.
*   Each block is annotated with a provenance line containing its original `source`, `job\_id`, `run\_id`, and original sample index (`orig\_s`) to prevent data loss or duplicate counting.
*   Stale files from previous runs are cleaned up automatically.

### Stage 4: Process and Analyze
Run the central post-processing bash wrapper:

    ./python/run\_pipeline.sh

This script runs the following steps sequentially [1]:
1.  Executes `merge\_samples.py` (Stage 3).
2.  Runs `process\_data.py` to create a master database (`long\_format.csv`) containing every individual sample and its provenance, alongside pivoted summary tables (`pivot\_mean\_p.csv`, `pivot\_mean\_t.csv`, `pivot\_sample\_count.csv`, `pivot\_solved\_rate.csv`).
3.  Runs `analyze\_data.py` to compute pairwise success probability ratios (`ratio.csv`) and find the individual samples producing the largest absolute separation between any two protocols (`max\_separation.csv`).

All resulting files are written directly to `results/processed/unitary/` or `results/processed/channels/`.

### Stage 5: Visualize
The repository includes a zero-dependency, browser-based plotting tool, `channel\_bench.html` (built with vanilla JS and Chart.js).
1.  Open `channel\_bench.html` in any web browser.
2.  Drag and drop any of the processed pivot CSV files (e.g., `pivot\_mean\_p.csv`) or `long\_format.csv` directly into the window.
3.  Use the interface to choose the independent variable, toggle specific protocol traces, select logarithmic scaling, apply line-smoothing, or adjust aesthetics for colorblind-friendly academic publishing.
4.  Export vector graphics directly as `.svg` files.

---

## Data Traceability

Every final result in `results/processed/` is fully auditable [1]. The pipeline preserves data lineage at every step [1]:

`Pivot Table / Analysis File` $\rightarrow$ `long\_format.csv` $\rightarrow$ `Data/<class>/<N>\_<n>.txt` $\rightarrow$ `slurm.out` $\rightarrow$ `results/raw/<RUN\_ID>/`

By looking up the provenance comment of any sample inside `Data/`, you can retrieve the exact random seed, physical parameters, and generated Choi matrices (`channels\_used\_diff.mat`) that produced the data point.

---

## Dependencies

*   **MATLAB** (R2020a or later recommended)
*   **CVX** (configured with the **Mosek** solver)
*   **QETLAB** (for `PartialTrace` and `PermuteSystems` functions)
*   **Python 3.x** (with `numpy` and `pandas` libraries)
