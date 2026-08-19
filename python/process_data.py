#!/usr/bin/env python3
"""
process_data.py

Reads the consolidated .txt files that merge_samples.py writes into
Data/unitary/ and Data/channels/ (one file per N, e.g. "3_10.txt"),
and turns them into:

  1. A single long-format table with one row per (N, sample, protocol),
     carrying full provenance (source .out file, job id, run id) so any
     number can be traced back to results/raw/<run_id>/.

  2. Pivoted N-by-protocol summary tables, written to results/processed/,
     computed separately for the 'unitary' and 'channels' data:
       - pivot_mean_p.csv      average success probability
       - pivot_mean_t.csv      average solve time (seconds), plus a
                                 total_time column (row-wise sum)
       - pivot_sample_count.csv  how many samples contributed to each cell
       - pivot_solved_rate.csv  fraction of samples with status "Solved"

Usage:
    python3 process_data.py
    python3 process_data.py --data-dir /path/to/Data --output-dir /path/to/results/processed
"""

import argparse
import csv
import glob
import os
import re
import sys
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
DEFAULT_DATA_DIR = os.path.join(REPO_ROOT, 'Data')
DEFAULT_OUTPUT_DIR = os.path.join(REPO_ROOT, 'results', 'processed')

# "N= 3 s=1/10 [unitary|diff]"
HEADER_RE = re.compile(r'^N=\s*(\d+)\s+s=\s*(\d+)/(\d+)\s+\[([^|\]]+)\|([^\]]+)\]')

# "  # source=unitary_upper_bound_N3_12178092.out job_id=12178092 run_id=... orig_s=1/5"
# (optional -- older merged files predate this line)
PROV_RE = re.compile(
    r'^\s*#\s*source=(\S+)\s+job_id=(\S+)\s+run_id=(\S+)\s+orig_s=(\d+)/(\d+)'
)

# "  PAR        p=1.000000 t=7.2s [Solved]"
RESULT_RE = re.compile(r'^\s*(\S+)\s+p=([\-\d.]+)\s+t=([\-\d.]+)s\s+\[(\w+)\]')


def parse_merged_file(filepath):
    """Yields one dict per (sample, protocol) result row."""
    with open(filepath, 'r') as f:
        lines = f.readlines()

    rows = []
    cur = None  # current block's shared metadata

    for line in lines:
        stripped = line.rstrip('\n')

        m = HEADER_RE.match(stripped)
        if m:
            cur = {
                'N': int(m.group(1)),
                'sample': int(m.group(2)),
                'n_samples_in_block': int(m.group(3)),
                'channel_type': m.group(4).strip(),
                'mode': m.group(5).strip(),
                'source_file': 'unknown',
                'job_id': 'unknown',
                'run_id': 'unknown',
                'orig_s': '',
                'orig_total': '',
            }
            continue

        pm = PROV_RE.match(stripped)
        if pm and cur is not None:
            cur['source_file'] = pm.group(1)
            cur['job_id'] = pm.group(2)
            cur['run_id'] = pm.group(3)
            cur['orig_s'] = pm.group(4)
            cur['orig_total'] = pm.group(5)
            continue

        rm = RESULT_RE.match(stripped)
        if rm and cur is not None:
            protocol, p, t, status = rm.groups()
            rows.append(dict(
                cur,
                protocol=protocol,
                p=float(p),
                t=float(t),
                status=status,
            ))

    return rows


def collect_long_format(data_dir, channel_class):
    """channel_class is 'unitary' or 'channels' -- the Data/<subfolder> name."""
    folder = os.path.join(data_dir, channel_class)
    if not os.path.isdir(folder):
        return []

    all_rows = []
    for path in sorted(glob.glob(os.path.join(folder, '*.txt'))):
        rows = parse_merged_file(path)
        for r in rows:
            r['channel_class'] = channel_class
            r['merged_file'] = os.path.basename(path)
        all_rows.extend(rows)
    return all_rows


LONG_FORMAT_COLUMNS = [
    'channel_class', 'merged_file', 'source_file', 'job_id', 'run_id',
    'channel_type', 'mode', 'N', 'sample', 'n_samples_in_block',
    'orig_s', 'orig_total', 'protocol', 'p', 't', 'status',
]


def write_long_format(rows, out_path):
    with open(out_path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=LONG_FORMAT_COLUMNS)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r.get(k, '') for k in LONG_FORMAT_COLUMNS})


def ordered_protocols(rows):
    seen = []
    for r in rows:
        if r['protocol'] not in seen:
            seen.append(r['protocol'])
    return seen


def build_pivots(rows):
    """
    Returns (protocols, pivot_mean_p, pivot_mean_t, pivot_sample_count,
    pivot_solved_rate), each of the pivot_* dicts keyed by N -> {protocol: value}.
    """
    protocols = ordered_protocols(rows)

    by_np = defaultdict(list)  # (N, protocol) -> list of rows
    for r in rows:
        by_np[(r['N'], r['protocol'])].append(r)

    Ns = sorted(set(r['N'] for r in rows))

    pivot_mean_p = {N: {} for N in Ns}
    pivot_mean_t = {N: {} for N in Ns}
    pivot_sample_count = {N: {} for N in Ns}
    pivot_solved_rate = {N: {} for N in Ns}

    for N in Ns:
        for proto in protocols:
            group = by_np.get((N, proto), [])
            if not group:
                continue
            n = len(group)
            pivot_mean_p[N][proto] = sum(g['p'] for g in group) / n
            pivot_mean_t[N][proto] = sum(g['t'] for g in group) / n
            pivot_sample_count[N][proto] = n
            solved = sum(1 for g in group if g['status'] == 'Solved')
            pivot_solved_rate[N][proto] = solved / n

    return protocols, pivot_mean_p, pivot_mean_t, pivot_sample_count, pivot_solved_rate


def write_pivot_csv(protocols, pivot, out_path, fmt='{:.6f}', extra_total_time=False):
    header = ['N'] + protocols
    if extra_total_time:
        header = header + ['total_time']
    with open(out_path, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        for N in sorted(pivot.keys()):
            vals = pivot[N]
            row = [N]
            row_vals = []
            for proto in protocols:
                v = vals.get(proto)
                row_vals.append(v)
                row.append(fmt.format(v) if v is not None else '')
            if extra_total_time:
                total = sum(v for v in row_vals if v is not None)
                row.append(fmt.format(total))
            writer.writerow(row)


def process_channel_class(data_dir, output_dir, channel_class):
    rows = collect_long_format(data_dir, channel_class)
    if not rows:
        print(f"  [{channel_class}] no data found under {os.path.join(data_dir, channel_class)}, skipping.")
        return

    class_dir = os.path.join(output_dir, channel_class)
    os.makedirs(class_dir, exist_ok=True)

    long_path = os.path.join(class_dir, 'long_format.csv')
    write_long_format(rows, long_path)
    print(f"  [{channel_class}] wrote {long_path} ({len(rows)} rows)")

    protocols, mean_p, mean_t, sample_count, solved_rate = build_pivots(rows)

    write_pivot_csv(protocols, mean_p, os.path.join(class_dir, 'pivot_mean_p.csv'))
    write_pivot_csv(protocols, mean_t, os.path.join(class_dir, 'pivot_mean_t.csv'), extra_total_time=True)
    write_pivot_csv(protocols, sample_count, os.path.join(class_dir, 'pivot_sample_count.csv'), fmt='{:.0f}')
    write_pivot_csv(protocols, solved_rate, os.path.join(class_dir, 'pivot_solved_rate.csv'))
    print(f"  [{channel_class}] wrote pivot_mean_p.csv, pivot_mean_t.csv, "
          f"pivot_sample_count.csv, pivot_solved_rate.csv -> {class_dir}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', default=DEFAULT_DATA_DIR,
                         help="Folder containing 'unitary' and 'channels' subfolders (default: Data/)")
    parser.add_argument('--output-dir', default=DEFAULT_OUTPUT_DIR,
                         help="Where to write long_format.csv and pivot_*.csv (default: results/processed)")
    args = parser.parse_args()

    if not os.path.isdir(args.data_dir):
        sys.exit(f"Error: data directory not found: {args.data_dir}")

    os.makedirs(args.output_dir, exist_ok=True)

    print(f"Processing {args.data_dir} -> {args.output_dir}")
    for channel_class in ('unitary', 'channels'):
        process_channel_class(args.data_dir, args.output_dir, channel_class)

    print("Done.")


if __name__ == '__main__':
    main()