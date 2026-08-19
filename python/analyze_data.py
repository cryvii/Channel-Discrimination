#!/usr/bin/env python3
"""
analyze_data.py

Reads the output of process_data.py (results/processed/<channel_class>/) and
computes two further summaries, for each of 'unitary' and 'channels':

  1. ratio.csv        For every N, the ratio of every protocol's mean
                       probability to every other protocol's (both
                       directions), from pivot_mean_p.csv.

  2. max_separation.csv
                       For every N and every pair of protocols, the single
                       sample where the gap between them is largest --
                       including that sample's run_id/job_id so you can
                       go straight to results/raw/<run_id>/ and see the
                       exact channels that produced it.

Usage:
    python3 analyze_data.py
    python3 analyze_data.py --processed-dir /path/to/results/processed
"""

import argparse
import csv
import itertools
import os
import sys
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
DEFAULT_PROCESSED_DIR = os.path.join(REPO_ROOT, 'results', 'processed')


def read_pivot_mean_p(path):
    """Returns (protocols, {N: {protocol: value}})."""
    with open(path, 'r', newline='') as fh:
        reader = csv.reader(fh)
        header = next(reader, None)
        if not header or header[0] != 'N':
            return [], {}
        protocols = header[1:]
        data = {}
        for row in reader:
            if not row:
                continue
            N = int(row[0])
            vals = {}
            for proto, v in zip(protocols, row[1:]):
                if v != '':
                    vals[proto] = float(v)
            data[N] = vals
    return protocols, data


def write_ratio_csv(protocols, data, out_path):
    pairs = list(itertools.permutations(protocols, 2))  # (denom, numer)
    header = ['N'] + ['{}/{}'.format(numer, denom) for denom, numer in pairs]
    with open(out_path, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        for N in sorted(data.keys()):
            vals = data[N]
            row = [N]
            for denom, numer in pairs:
                d = vals.get(denom)
                n = vals.get(numer)
                if d is None or n is None or d == 0:
                    row.append('')
                else:
                    row.append('{:.6f}'.format(n / d))
            writer.writerow(row)


def read_long_format(path):
    with open(path, 'r', newline='') as fh:
        reader = csv.DictReader(fh)
        return list(reader)


def find_max_separations(protocols, rows):
    """
    For every (lower, higher) pair of protocols (in the given order),
    for every N, find the single sample with the largest (higher - lower)
    gap in success probability. Returns {N: {(lower, higher): row_dict}}
    where row_dict has 'diff', plus the winning sample's provenance.
    """
    # index: (N, sample, protocol) -> row
    by_key = {}
    samples_by_N = defaultdict(set)
    for r in rows:
        N = int(r['N'])
        sample = int(r['sample'])
        by_key[(N, sample, r['protocol'])] = r
        samples_by_N[N].add(sample)

    pairs = list(itertools.combinations(protocols, 2))  # (lower, higher)
    result = defaultdict(dict)

    for N, samples in samples_by_N.items():
        for lower, higher in pairs:
            best = None
            for sample in samples:
                lo = by_key.get((N, sample, lower))
                hi = by_key.get((N, sample, higher))
                if lo is None or hi is None:
                    continue
                diff = float(hi['p']) - float(lo['p'])
                if best is None or diff > best['diff']:
                    best = {
                        'diff': diff,
                        'sample': sample,
                        'run_id': hi.get('run_id', 'unknown'),
                        'job_id': hi.get('job_id', 'unknown'),
                        'source_file': hi.get('source_file', 'unknown'),
                        'orig_s': hi.get('orig_s', ''),
                    }
            if best is not None:
                result[N][(lower, higher)] = best

    return pairs, result


def write_max_separation_csv(pairs, n_to_best, out_path):
    header = ['N']
    for lower, higher in pairs:
        label = '{}-{}'.format(higher, lower)
        header += [label + '_max_diff', label + '_sample',
                   label + '_run_id', label + '_job_id']
    with open(out_path, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        for N in sorted(n_to_best.keys()):
            best = n_to_best[N]
            row = [N]
            for pair in pairs:
                b = best.get(pair)
                if b is None:
                    row += ['', '', '', '']
                else:
                    row += ['{:.6f}'.format(b['diff']), b['sample'],
                            b['run_id'], b['job_id']]
            writer.writerow(row)


def process_channel_class(processed_dir, channel_class):
    class_dir = os.path.join(processed_dir, channel_class)
    pivot_path = os.path.join(class_dir, 'pivot_mean_p.csv')
    long_path = os.path.join(class_dir, 'long_format.csv')

    if not os.path.isfile(pivot_path) or not os.path.isfile(long_path):
        print(f"  [{channel_class}] no processed data found in {class_dir}, skipping "
              f"(run process_data.py first).")
        return

    protocols, pivot_data = read_pivot_mean_p(pivot_path)
    if protocols and pivot_data:
        out_path = os.path.join(class_dir, 'ratio.csv')
        write_ratio_csv(protocols, pivot_data, out_path)
        print(f"  [{channel_class}] wrote {out_path}")

    rows = read_long_format(long_path)
    if rows:
        pairs, n_to_best = find_max_separations(protocols, rows)
        out_path = os.path.join(class_dir, 'max_separation.csv')
        write_max_separation_csv(pairs, n_to_best, out_path)
        print(f"  [{channel_class}] wrote {out_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--processed-dir', default=DEFAULT_PROCESSED_DIR,
                         help="Output folder from process_data.py (default: results/processed)")
    args = parser.parse_args()

    if not os.path.isdir(args.processed_dir):
        sys.exit(f"Error: processed directory not found: {args.processed_dir}\n"
                  f"Run process_data.py first.")

    print(f"Analyzing {args.processed_dir}")
    for channel_class in ('unitary', 'channels'):
        process_channel_class(args.processed_dir, channel_class)

    print("Done.")


if __name__ == '__main__':
    main()