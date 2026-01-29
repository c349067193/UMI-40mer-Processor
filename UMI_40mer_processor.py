#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Spill-to-disk combiner for 20+20 nt tags with QUALITY SPLIT.
- Combines R1[:20] + {revcomp|comp|raw}(R2[:20]) if both reads match the 30-nt motifs at pos 21–50.
- Splits counts into two outputs:
    * Q-pass (all 20+20 bases have Phred >= --q-threshold, default 30)
    * Q-fail (anything else)
- Flushes partial counts to disk and externally merges (exact unique totals, bounded RAM).
- Multi-lane globbing and optional pigz decompression.

Example (exact positions, your motifs):
  python3 UMI_40mer_processor.py 'BMR02b_L*_1.fq.gz' \
    -o BMR02b_20nt_exact \
    --fwd-motif ggcattctggatagtgtcaagggaataaat \
    --rev-motif gcaaatatcatttattcccttgacactatc \
    --r2-head revcomp \
    --r1-scan 0 --r2-scan 0 \
    --q-threshold 30 \
    --jobs 4 --pigz-threads 1 \
    --progress-every 1000000 \
    --spill-dir BMR02b_spill \
    --flush-every-unique 2000000
"""

from __future__ import annotations

import argparse, gzip, io, os, re, sys, glob, shutil, subprocess, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Iterable, Tuple, List, Dict

# -------------------------- IO helpers --------------------------

def is_gz(path: str) -> bool:
    return path.endswith(".gz")

def pigz_available() -> bool:
    return shutil.which("pigz") is not None

def open_text(path: str, pigz_threads: int = 0):
    """Open text file or .gz (optionally via pigz)."""
    if is_gz(path) and pigz_threads > 0 and pigz_available():
        proc = subprocess.Popen(["pigz", "-dc", "-p", str(pigz_threads), path], stdout=subprocess.PIPE)
        return io.TextIOWrapper(proc.stdout, encoding="utf-8", newline="")
    if is_gz(path):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", newline="")
    return open(path, "rt", encoding="utf-8", newline="")

def fastq_read_iter(path: str, pigz_threads: int = 0) -> Iterable[Tuple[str, str]]:
    """Yield (seq, qual) from FASTQ(.gz)."""
    with open_text(path, pigz_threads=pigz_threads) as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            seq = fh.readline().rstrip("\n")
            plus = fh.readline()
            qual = fh.readline().rstrip("\n")
            if not qual:
                break
            yield seq, qual

# -------------------------- Bio utils --------------------------

def revcomp(seq: str) -> str:
    tbl = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(tbl)[::-1]

def comp(seq: str) -> str:
    tbl = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(tbl)

def matches_motif_with_offset(seq: str, motif: str, center_pos_1based: int = 21, scan: int = 0):
    """Case-insensitive; return (ok, start_pos_1based)."""
    if not motif:
        return False, -1
    s = seq.upper()
    m = motif.upper()
    L = len(m)
    start0 = center_pos_1based - 1
    lo = max(0, start0 - scan)
    hi = min(start0 + scan, len(s) - L)
    if hi < lo:
        return False, -1
    for off in range(lo, hi + 1):
        if s[off : off + L] == m:
            return True, off + 1
    return False, -1

def all_q_ge(qual: str, n: int, thr: int) -> bool:
    """Check first n bases have Phred >= thr (Sanger/Illumina+33)."""
    # Avoid slow Python loops on long strings; but here n=20, so a loop is fine/readable.
    for c in qual[:n]:
        if (ord(c) - 33) < thr:
            return False
    return True

# --------------------- Spill/flush + merge ---------------------

def flush_counts_chunk(counts: Dict[str, int], lane_tag: str, cls: str, spill_dir: str, chunk_idx: int) -> str:
    """Write one chunk file for a class (Qpass / Qfail)."""
    if not counts:
        return ""
    os.makedirs(spill_dir, exist_ok=True)
    # cls is 'Q{thr}' or 'subQ{thr}'
    path = os.path.join(spill_dir, f"{lane_tag}.{cls}.chunk{chunk_idx:05d}.tsv")
    with open(path, "w", encoding="utf-8") as out:
        for seq, c in counts.items():
            out.write(f"{seq}\t{c}\n")
    return path

def external_merge_counts(chunk_paths: List[str], out_counts_path: str, temp_dir: str) -> int:
    """Sort+merge (seq\tcount) chunks into a final file. Returns exact #unique sequences written."""
    unique_written = 0
    if not chunk_paths:
        open(out_counts_path, "w").close()
        return unique_written

    sort_bin = shutil.which("gsort") or "sort"
    sort_cmd = [sort_bin]
    # If GNU sort available, consider parallel + memory; tune if desired:
    if os.path.basename(sort_bin) == "gsort":
        # Conservatively use 4 threads and 4G; adjust to your machine.
        sort_cmd += ["--parallel", "4", "-S", "4G"]
    if temp_dir:
        sort_cmd += ["-T", temp_dir]

    env = os.environ.copy()
    env["LC_ALL"] = "C"

    sort_proc = subprocess.Popen(sort_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, env=env, text=True, bufsize=1)

    # Stream all chunks into sort's stdin
    for p in chunk_paths:
        with open(p, "r", encoding="utf-8") as fh:
            for line in fh:
                sort_proc.stdin.write(line)
    sort_proc.stdin.close()

    # Stream-aggregate sorted lines
    with open(out_counts_path, "w", encoding="utf-8") as out:
        prev_key = None
        total = 0
        for line in sort_proc.stdout:
            line = line.rstrip("\n")
            if not line:
                continue
            try:
                key, val_s = line.split("\t", 1)
                val = int(val_s)
            except ValueError:
                continue
            if prev_key is None:
                prev_key = key
                total = val
            elif key == prev_key:
                total += val
            else:
                out.write(f"{prev_key}\t{total}\n")
                unique_written += 1
                prev_key = key
                total = val
        if prev_key is not None:
            out.write(f"{prev_key}\t{total}\n")
            unique_written += 1

    sort_proc.wait()
    return unique_written

# ----------------------- Pair processing -----------------------

def process_one_pair(
    r1: str,
    r2: str,
    fwd_motif: str,
    rev_motif: str,
    r2_head_mode: str,
    pigz_threads: int,
    progress_every: int,
    r1_scan: int,
    r2_scan: int,
    max_pairs_per_lane: int,
    flush_every_unique: int,
    spill_dir: str,
    q_threshold: int,
) -> Tuple[List[str], List[str], int, int, int, int, int, int, int]:
    """
    Process a single R1/R2 pair.
    Returns:
      (chunk_paths_pass, chunk_paths_fail,
       total_pairs, matched_pairs, matched_pass, matched_fail,
       nonmatch_pairs, too_short, r1_fail, r2_fail)
    """
    total_pairs = matched_pairs = matched_pass = matched_fail = 0
    nonmatch_pairs = too_short = r1_fail = r2_fail = 0

    rc = revcomp
    cp = comp
    head_mode = r2_head_mode

    it1 = fastq_read_iter(r1, pigz_threads)
    it2 = fastq_read_iter(r2, pigz_threads)

    t0 = time.time()
    lane_tag = os.path.basename(r1)

    counts_pass: Dict[str, int] = {}
    counts_fail: Dict[str, int] = {}
    chunk_idx_pass = chunk_idx_fail = 0
    chunks_pass: List[str] = []
    chunks_fail: List[str] = []

    for (s1, q1), (s2, q2) in zip(it1, it2):
        total_pairs += 1
        if len(s1) < 50 or len(s2) < 50:
            nonmatch_pairs += 1
            too_short += 1
        else:
            ok1, _ = matches_motif_with_offset(s1, fwd_motif, center_pos_1based=21, scan=r1_scan)
            ok2, _ = matches_motif_with_offset(s2, rev_motif, center_pos_1based=21, scan=r2_scan)
            if ok1 and ok2:
                matched_pairs += 1
                head2 = s2[:20]
                if head_mode == "revcomp":
                    head2 = rc(head2)
                elif head_mode == "comp":
                    head2 = cp(head2)
                combined = s1[:20] + head2

                # QC on R1[:20] and R2[:20] qualities
                pass_q = all_q_ge(q1, 20, q_threshold) and all_q_ge(q2, 20, q_threshold)
                if pass_q:
                    counts_pass[combined] = counts_pass.get(combined, 0) + 1
                    matched_pass += 1
                    if flush_every_unique and len(counts_pass) >= flush_every_unique:
                        p = flush_counts_chunk(counts_pass, lane_tag, f"Q{q_threshold}", spill_dir, chunk_idx_pass)
                        if p:
                            chunks_pass.append(p)
                        counts_pass.clear()
                        chunk_idx_pass += 1
                else:
                    counts_fail[combined] = counts_fail.get(combined, 0) + 1
                    matched_fail += 1
                    if flush_every_unique and len(counts_fail) >= flush_every_unique:
                        p = flush_counts_chunk(counts_fail, lane_tag, f"subQ{q_threshold}", spill_dir, chunk_idx_fail)
                        if p:
                            chunks_fail.append(p)
                        counts_fail.clear()
                        chunk_idx_fail += 1
            else:
                nonmatch_pairs += 1
                if not ok1:
                    r1_fail += 1
                if not ok2:
                    r2_fail += 1

        if max_pairs_per_lane and total_pairs >= max_pairs_per_lane:
            break

        if progress_every > 0 and (total_pairs % progress_every == 0):
            elapsed = max(time.time() - t0, 1e-9)
            rate = total_pairs / elapsed
            print(
                f"[PROGRESS] {lane_tag}: {total_pairs:,} pairs | matched {matched_pairs:,} "
                f"(Q≥{q_threshold}: {matched_pass:,} / subQ: {matched_fail:,}) | non-match {nonmatch_pairs:,} | ~{rate:,.0f} pairs/s",
                file=sys.stderr, flush=True
            )

    # final flush (both classes)
    if counts_pass:
        p = flush_counts_chunk(counts_pass, lane_tag, f"Q{q_threshold}", spill_dir, chunk_idx_pass)
        if p: chunks_pass.append(p)
        counts_pass.clear()
    if counts_fail:
        p = flush_counts_chunk(counts_fail, lane_tag, f"subQ{q_threshold}", spill_dir, chunk_idx_fail)
        if p: chunks_fail.append(p)
        counts_fail.clear()

    elapsed = max(time.time() - t0, 1e-9)
    rate = total_pairs / elapsed
    print(
        f"[PROGRESS] {lane_tag}: {total_pairs:,} pairs processed (final) | matched {matched_pairs:,} "
        f"(Q≥{q_threshold}: {matched_pass:,} / subQ: {matched_fail:,}) | non-match {nonmatch_pairs:,} | ~{rate:,.0f} pairs/s",
        file=sys.stderr, flush=True
    )

    return chunks_pass, chunks_fail, total_pairs, matched_pairs, matched_pass, matched_fail, nonmatch_pairs, too_short, r1_fail, r2_fail

# ----------------------- File pair helpers ---------------------

_R1_R2_REPLACEMENTS = [
    (re.compile(r'([._-])1(\.(?:f(?:ast)?q)(?:\.gz)?)$', re.IGNORECASE), r'\1 2\2'),
    (re.compile(r'([/_-])R1([/_-])', re.IGNORECASE), r'\1R2\2'),
    (re.compile(r'([._-])R1(\.(?:f(?:ast)?q)(?:\.gz)?)$', re.IGNORECASE), r'\1R2\2'),
    (re.compile(r'(?:^|[/_\.-])L(\d+)_1(\.(?:f(?:ast)?q)(?:\.gz)?)$', re.IGNORECASE), r'L\1_2\2'),
]

def guess_r2_path(r1_path: str) -> str | None:
    dirname, fname = os.path.split(r1_path)
    cands = []
    for pat, repl in _R1_R2_REPLACEMENTS:
        new = pat.sub(repl, fname)
        if new != fname:
            cands.append(os.path.join(dirname, new.replace(' 2', '2')))
    for a, b in [('_1.fq.gz','_2.fq.gz'),('_1.fastq.gz','_2.fastq.gz'),
                 ('_R1_','_R2_'),('_R1.fq.gz','_R2.fq.gz'),('_R1.fastq.gz','_R2.fastq.gz'),
                 ('_1.fq','_2.fq'),('_1.fastq','_2.fastq')]:
        if fname.endswith(a):
            cands.append(os.path.join(dirname, fname[:-len(a)] + b))
    for cand in cands:
        if os.path.exists(cand):
            return cand
    return None

def expand_r1_list(patterns: List[str]) -> List[str]:
    r1_files: List[str] = []
    for pat in patterns:
        expanded = sorted(glob.glob(pat)) if any(ch in pat for ch in "*?[") else [pat]
        r1_files.extend(expanded)
    seen = set(); uniq: List[str] = []
    for p in r1_files:
        if p not in seen:
            uniq.append(p); seen.add(p)
    return uniq

# -------------------------- Orchestration ----------------------

def process_pairs_spill_qc(
    r1_files: List[str],
    fwd_motif: str,
    rev_motif: str,
    r2_head_mode: str,
    out_prefix: str,
    jobs: int,
    pigz_threads: int,
    progress_every: int,
    r1_scan: int,
    r2_scan: int,
    max_pairs_per_lane: int,
    flush_every_unique: int,
    spill_dir: str,
    q_threshold: int,
) -> None:
    pairs: List[Tuple[str, str]] = []
    for r1 in r1_files:
        r2 = guess_r2_path(r1)
        if r2 is None:
            print(f"[WARN] Could not locate matching R2 for {r1}; skipping.", file=sys.stderr)
            continue
        pairs.append((r1, r2))

    if not pairs:
        print("[ERROR] No valid R1/R2 pairs found.", file=sys.stderr)
        sys.exit(2)

    # workers
    try:
        import multiprocessing as mp
        cpu_cnt = mp.cpu_count()
    except Exception:
        cpu_cnt = 1
    if jobs is None or jobs <= 0:
        jobs = min(len(pairs), cpu_cnt)

    print(f"[INFO] Processing {len(pairs)} pair(s) with {jobs} worker(s).", file=sys.stderr)
    if pigz_threads > 0:
        if not pigz_available():
            print("[WARN] pigz not found; falling back to Python gzip. Install via: brew install pigz", file=sys.stderr)
        else:
            print(f"[INFO] pigz decompression enabled with {pigz_threads} thread(s) per worker.", file=sys.stderr)

    os.makedirs(spill_dir, exist_ok=True)

    all_chunks_pass: List[str] = []
    all_chunks_fail: List[str] = []
    total = matched = matched_pass = matched_fail = nonmatch = 0
    too_short = r1_fail = r2_fail = 0

    with ProcessPoolExecutor(max_workers=jobs) as ex:
        futs = [ex.submit(
            process_one_pair,
            r1, r2,
            fwd_motif, rev_motif, r2_head_mode,
            pigz_threads, progress_every, r1_scan, r2_scan, max_pairs_per_lane,
            flush_every_unique, spill_dir, q_threshold
        ) for (r1, r2) in pairs]

        for fut in as_completed(futs):
            chunks_p, chunks_f, t, m, mp_, mf_, n, ts, r1f, r2f = fut.result()
            all_chunks_pass.extend(chunks_p)
            all_chunks_fail.extend(chunks_f)
            total += t; matched += m; matched_pass += mp_; matched_fail += mf_
            nonmatch += n; too_short += ts; r1_fail += r1f; r2_fail += r2f

    # Merge both classes
    counts_pass_path = f"{out_prefix}_combined_counts.Q{q_threshold}.txt"
    counts_fail_path = f"{out_prefix}_combined_counts.subQ{q_threshold}.txt"
    print(f"[INFO] Merging {len(all_chunks_pass)} Q≥{q_threshold} chunk(s) into {counts_pass_path} ...", file=sys.stderr)
    uniq_pass = external_merge_counts(all_chunks_pass, counts_pass_path, temp_dir=spill_dir)
    print(f"[INFO] Merging {len(all_chunks_fail)} subQ chunk(s) into {counts_fail_path} ...", file=sys.stderr)
    uniq_fail = external_merge_counts(all_chunks_fail, counts_fail_path, temp_dir=spill_dir)

    # Write summary
    summary_path = f"{out_prefix}_summary.txt"
    with open(summary_path, "w", encoding="utf-8") as out:
        out.write(f"Total paired reads: {total}\n")
        out.write(f"Combined (20+20) paired reads: {matched}\n")
        out.write(f"  Q≥{q_threshold} pairs: {matched_pass}\n")
        out.write(f"  subQ pairs: {matched_fail}\n")
        out.write(f"Non-matching pairs: {nonmatch}\n")
        out.write("Non-matching breakdown:\n")
        out.write(f"  too_short(<50nt): {too_short}\n")
        out.write(f"  R1 motif missing (±{r1_scan}): {r1_fail}\n")
        out.write(f"  R2 motif missing (±{r2_scan}): {r2_fail}\n")
        out.write("Unique 40nt sequences (exact):\n")
        out.write(f"  Q≥{q_threshold}: {uniq_pass}\n")
        out.write(f"  subQ: {uniq_fail}\n")
        out.write("Parameters:\n")
        out.write(f"  r2_head_mode: {r2_head_mode}\n")
        out.write(f"  fwd_motif: {fwd_motif}\n")
        out.write(f"  rev_motif: {rev_motif}\n")
        out.write(f"  r1_scan: ±{r1_scan}\n")
        out.write(f"  r2_scan: ±{r2_scan}\n")
        out.write(f"  pigz_threads: {pigz_threads}\n")
        out.write(f"  progress_every: {progress_every}\n")
        out.write(f"  flush_every_unique: {flush_every_unique}\n")
        out.write(f"  spill_dir: {spill_dir}\n")
        out.write(f"  q_threshold: {q_threshold}\n")

    print(f"[DONE] Wrote:\n  {counts_pass_path}\n  {counts_fail_path}\n  {summary_path}", file=sys.stderr)

# ------------------------------- CLI -------------------------------

def main():
    p = argparse.ArgumentParser(description="Combine 20+20 nt tags with quality split (spill-to-disk).")
    p.add_argument("r1_files", nargs="+", help='R1 files or globs (e.g., "BMR*_L*_1.fq.gz")')
    p.add_argument("-o", "--out-prefix", default="combined_20nt", help="Output prefix")
    p.add_argument("--fwd-motif", required=True, help="30nt motif expected near R1[21:50]")
    p.add_argument("--rev-motif", required=True, help="30nt motif expected near R2[21:50] (RAW R2 orientation)")
    p.add_argument("--r2-head", dest="r2_head_mode", choices=["revcomp", "comp", "raw"], default="revcomp",
                   help="Transform R2[:20] before concatenation (default: revcomp)")
    p.add_argument("--jobs", type=int, default=0, help="Worker processes (default: min(#pairs, CPU cores))")
    p.add_argument("--pigz-threads", type=int, default=0, help="If >0 and .gz, use pigz with N threads")
    p.add_argument("--progress-every", type=int, default=1_000_000, help="Progress interval per lane (0=off)")
    p.add_argument("--r1-scan", type=int, default=0, help="Allow ±N nt shift for R1 motif (default 0)")
    p.add_argument("--r2-scan", type=int, default=0, help="Allow ±N nt shift for R2 motif (default 0)")
    p.add_argument("--max-pairs-per-lane", type=int, default=0, help="Cap per-lane reads for testing (0=all)")
    p.add_argument("--flush-every-unique", type=int, default=2_000_000, help="Flush counts when UNIQUE hits N")
    p.add_argument("--spill-dir", default=None, help="Spill directory (default: <out-prefix>_spill)")
    p.add_argument("--q-threshold", type=int, default=30, help="Quality threshold for Q-pass (default 30)")

    args = p.parse_args()

    r1_list = expand_r1_list(args.r1_files)
    if not r1_list:
        print("[ERROR] No R1 files matched.", file=sys.stderr)
        sys.exit(2)

    spill_dir = args.spill_dir or f"{args.out_prefix}_spill"

    if len(args.fwd_motif) != 30 or len(args.rev_motif) != 30:
        print(f"[WARN] Expected 30-nt motifs; got {len(args.fwd_motif)} and {len(args.rev_motif)}", file=sys.stderr)

    process_pairs_spill_qc(
        r1_list,
        args.fwd_motif,
        args.rev_motif,
        args.r2_head_mode,
        args.out_prefix,
        args.jobs,
        args.pigz_threads,
        args.progress_every,
        args.r1_scan,
        args.r2_scan,
        args.max_pairs_per_lane,
        args.flush_every_unique,
        spill_dir,
        args.q_threshold,
    )

if __name__ == "__main__":
    main()
