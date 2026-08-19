#!/usr/bin/env python3
import os
import re
import sys
import glob
import argparse
from collections import defaultdict

# -------------------------------------------------------------------------
# DYNAMIC PATH RESOLUTION
# This resolves your workspace folders relative to the location of this script,
# meaning it will run perfectly from ANY terminal directory context.
# -------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
DEFAULT_LOG_DIR = os.path.join(REPO_ROOT, 'matlab', 'slurm_logs')
DEFAULT_OUTPUT_DIR = os.path.join(REPO_ROOT, 'Data')

# Regex to find sample block headers: e.g., "N= 2 s=1/50 [unitary|diff]"
HEADER_RE = re.compile(r'^N=\s*(\d+)\s+s=\s*(\d+)/(\d+)\s+\[([^\]]+)\]')

# Matches the filenames this script itself produces, e.g. "3_8.txt", "12_150.txt"
OUTPUT_FILENAME_RE = re.compile(r'^(\d+)_(\d+)\.txt$')

# Provenance info printed near the top of every slurm .out file by run_job.sh:
#   "  Job: unitary_upper_bound_N3  (ID: 12178092)"
#   "  RUN_ID: 20260818_145406_job12178092_unitary_upper_bound_N3"
JOB_ID_RE = re.compile(r'Job:\s*\S+\s*\(ID:\s*(\d+)\)')
RUN_ID_RE = re.compile(r'^\s*RUN_ID:\s*(\S+)')


def parse_file_provenance(lines):
    """
    Scans the top of a .out file for the job ID and RUN_ID printed by
    run_job.sh. Returns (job_id, run_id), each 'unknown' if not found.
    """
    job_id = 'unknown'
    run_id = 'unknown'
    # These lines always appear in the first ~10 lines of the file.
    for line in lines[:20]:
        if job_id == 'unknown':
            m = JOB_ID_RE.search(line)
            if m:
                job_id = m.group(1)
        if run_id == 'unknown':
            m = RUN_ID_RE.match(line)
            if m:
                run_id = m.group(1)
        if job_id != 'unknown' and run_id != 'unknown':
            break
    return job_id, run_id


def parse_out_file(filepath):
    """
    Parses an individual .out log file, extracting all sample blocks.
    Each block is tagged with provenance (source file, job id, run id,
    and its original sample index within that file) so it can be traced
    back to results/raw/<run_id>/ after merging.
    Returns a list of dictionaries, each representing a single parsed block.
    """
    blocks = []
    if not os.path.exists(filepath):
        print(f"Warning: File not found: {filepath}", file=sys.stderr)
        return blocks

    with open(filepath, 'r') as f:
        lines = f.readlines()

    source_file = os.path.basename(filepath)
    job_id, run_id = parse_file_provenance(lines)

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        match = HEADER_RE.match(line)
        if match:
            n_val = int(match.group(1))
            orig_s = int(match.group(2))
            orig_total = int(match.group(3))
            mode_str = match.group(4)
            protocol_lines = []
            i += 1

            # Consume lines belonging to this block until we hit a blank line,
            # another header, or end-of-run markers.
            while i < len(lines):
                next_line = lines[i]
                next_line_stripped = next_line.strip()
                if (HEADER_RE.match(next_line_stripped) or
                        next_line_stripped == "" or
                        "Saved results" in next_line_stripped or
                        "=== DONE ===" in next_line_stripped):
                    break
                protocol_lines.append(next_line.rstrip())
                i += 1

            blocks.append({
                'n': n_val,
                'mode': mode_str,
                'protocols': protocol_lines,
                'source_file': source_file,
                'job_id': job_id,
                'run_id': run_id,
                'orig_s': orig_s,
                'orig_total': orig_total,
            })
            continue
        i += 1
    return blocks


def clear_stale_outputs(target_dir, n_val, keep_filename):
    """
    Removes any previously-generated '<n_val>_<count>.txt' files in target_dir
    that don't match the file we're about to write. This is what prevents old
    runs (with a different sample count) from being left behind alongside the
    new one, e.g. '3_8.txt' lingering next to a fresh '3_10.txt'.
    """
    if not os.path.isdir(target_dir):
        return

    pattern = os.path.join(target_dir, f"{n_val}_*.txt")
    for existing_path in glob.glob(pattern):
        existing_name = os.path.basename(existing_path)
        m = OUTPUT_FILENAME_RE.match(existing_name)
        if not m:
            continue  # not one of our generated files; leave it alone
        existing_n = int(m.group(1))
        if existing_n == n_val and existing_name != keep_filename:
            os.remove(existing_path)
            print(f"  Removed stale file: {existing_path}")


def main():
    parser = argparse.ArgumentParser(description="Consolidate multiple .out log files, automatically grouping by N and channel type.")
    parser.add_argument('--log-dir', type=str, default=DEFAULT_LOG_DIR, help="Directory containing standard .out logs.")
    parser.add_argument('--output-dir', type=str, default=DEFAULT_OUTPUT_DIR, help="Base directory to save consolidated results (creates /unitary and /channels).")
    parser.add_argument('--files', type=str, nargs='+', help="Specific file names to merge (optional). If empty, auto-scans the log-dir.")

    args = parser.parse_args()

    # Resolve target files
    target_files = []
    if args.files:
        target_files = [os.path.join(args.log_dir, f) if not os.path.isabs(f) else f for f in args.files]
    else:
        # Automatically find all .out files in log-dir
        if os.path.exists(args.log_dir):
            all_files = os.listdir(args.log_dir)
            matching_files = [f for f in all_files if f.endswith('.out')]
            target_files = [os.path.join(args.log_dir, f) for f in sorted(matching_files)]
        else:
            print(f"Error: Log directory does not exist: {args.log_dir}", file=sys.stderr)
            sys.exit(1)

    if not target_files:
        print(f"No .out log files found in: {args.log_dir}")
        sys.exit(0)

    print(f"Scanning {len(target_files)} log file(s)...")

    # Extract all blocks from target files
    all_blocks = []
    for filepath in target_files:
        blocks = parse_out_file(filepath)
        all_blocks.extend(blocks)

    if len(all_blocks) == 0:
        print("No sample block data found inside the target .out files.", file=sys.stderr)
        sys.exit(1)

    # Group blocks dynamically by (N, channel_class)
    # Keys will be tuples: (n_val, 'unitary' or 'channels')
    grouped_data = defaultdict(list)

    for block in all_blocks:
        # Determine channel class
        mode_parts = block['mode'].split('|')
        channel_class = mode_parts[0].strip().lower()

        type_str = 'unitary' if channel_class == 'unitary' else 'channels'
        n_val = block['n']

        grouped_data[(n_val, type_str)].append(block)

    # Write grouped data to files
    print(f"\nProcessing complete. Writing files to: {args.output_dir}")

    # Sort by channel type, then by numerical N value
    sorted_keys = sorted(grouped_data.keys(), key=lambda x: (x[1], x[0]))

    for (n_val, type_str) in sorted_keys:
        blocks = grouped_data[(n_val, type_str)]
        total_samples = len(blocks)

        # Ensure correct folder directory structure
        target_dir = os.path.join(args.output_dir, type_str)
        os.makedirs(target_dir, exist_ok=True)

        out_filename = f"{n_val}_{total_samples}.txt"
        out_path = os.path.join(target_dir, out_filename)

        # Remove any stale file(s) left over from a previous run with a
        # different sample count for this same N, e.g. an old "3_8.txt"
        # before we write the current "3_10.txt".
        clear_stale_outputs(target_dir, n_val, out_filename)

        with open(out_path, 'w') as outfile:
            for idx, block in enumerate(blocks, 1):
                # Main header, same format as before (parseable by HEADER_RE).
                header = f"N= {n_val} s={idx}/{total_samples} [{block['mode']}]"
                outfile.write(header + "\n")
                # Provenance line: lets you trace this exact sample back to
                # its source .out file and the results/raw/<run_id>/ folder
                # containing the actual channels used (channels_used_*.txt/.mat).
                outfile.write(
                    f"  # source={block['source_file']} "
                    f"job_id={block['job_id']} "
                    f"run_id={block['run_id']} "
                    f"orig_s={block['orig_s']}/{block['orig_total']}\n"
                )
                for proto_line in block['protocols']:
                    outfile.write(proto_line + "\n")
                outfile.write("\n")

        print(f"  [{type_str.upper()}] N = {n_val:2d}: Consolidated {total_samples:3d} samples -> {type_str}/{out_filename}")


if __name__ == '__main__':
    main()