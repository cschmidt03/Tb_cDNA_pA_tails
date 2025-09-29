#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filter a BAM by detecting a 3' adapter inside soft-clips (read space),
optionally classify rejects by scanning the whole read, and (for kept reads)
extract a bounded window adjacent to the adapter.

------------------------------------------------------------
KEEP RULE (always)
------------------------------------------------------------
- Inspect ONLY the soft-clips (CIGAR 'S'):
  • 5′ soft-clip is searched for the REVERSE-COMPLEMENT of the adapter.
  • 3′ soft-clip is searched for the FORWARD adapter.
- A match is valid if it is:
  • exact substring, OR
  • ≤ mm_max mismatches (Hamming distance), OR
  • exactly one deletion in the adapter (window length = len(adapter)-1)
  (No combined mismatch+indel.)
- KEEP the read iff there is **exactly one** adapter hit across both soft-clips.

------------------------------------------------------------
REJECT CLASSIFICATION (optional, --fullread-classify)
------------------------------------------------------------
If enabled, rejected reads are additionally scanned across the FULL read
(both orientations, same matching rules) to bucket them as:
  • removed_no_adapter_anywhere
  • removed_adapter_multiple_anywhere
  • removed_adapter_present_not_in_softclip
  • removed_softclip_multiple_hits   (rare: >1 across soft-clips)
If not enabled (fast mode), report only:
  • removed_no_adapter_softclip
  • removed_adapter_multiple_softclip

------------------------------------------------------------
ADJACENT SEQUENCE EXTRACTION (kept reads only)
------------------------------------------------------------
If --extract-n > 0 and --extract-out is set, for each kept read:
  • 3′ adapter (FWD in 3′ soft-clip):
      output = sc3[ max(0, s - N) : e ]        # upstream ≤N + adapter, as-is
  • 5′ adapter (RC in 5′ soft-clip):
      tmp    = sc5[ s : min(len(sc5), e + N) ]  # adapter_RC + ≤N downstream
      output = revcomp(tmp)                     # reverse-complement to FWD
The slice is hard-clipped to the soft-clip bounds (if less than N available,
you get fewer). Results are written to TSV:
  read_id  side  hit_mode  adapter_start  adapter_end  extracted
Where (adapter_start, adapter_end) are coordinates within the corresponding
soft-clip string.

------------------------------------------------------------
COMMON OPTIONS
------------------------------------------------------------
--adapter CTGTAGGCACCATCAAT     # 3' adapter (forward)
--mm-max 2                      # allow up to 2 mismatches (default 1)
--no-del                        # disable the “1 deletion” allowance
--min-mapq 1                    # default removes MAPQ=0; set 0 to keep all
--fullread-classify             # classify rejects by scanning full read
--extract-n 10 --extract-out kept_adjacent.tsv

Example:
  python filter_adapter_categorized.py \
    -i input.bam -o filtered.bam \
    --adapter CTGTAGGCACCATCAAT \
    --mm-max 2 \
    --fullread-classify \
    --extract-n 10 --extract-out filtered.extracted.tsv
"""

import argparse
from collections import Counter
import sys
import pysam

ADAPTER_DEFAULT = "CTGTAGGCACCATCAAT"

# ---------------- utilities ----------------

def revcomp(s):
    tr = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return s.translate(tr)[::-1]

def ham(a, b):
    return sum(x != y for x, y in zip(a, b))

def get_softclips(read):
    """
    Return (softclip_5p_seq, softclip_3p_seq) sliced from read.query_sequence
    using leading/trailing CIGAR 'S' (soft-clip) operations.
    """
    q = read.query_sequence
    if not q:
        return None, None
    ct = read.cigartuples or []
    sc5 = q[:ct[0][1]] if ct and ct[0][0] == 4 else None
    sc3 = q[-ct[-1][1]:] if ct and ct[-1][0] == 4 else None
    return sc5, sc3

def dedup_hits(hits, slop=1):
    """
    Merge nearly-identical loci (start/end within 'slop'); prefer longer span.
    hits: list of (start, end, mode)
    """
    if not hits:
        return []
    hits = sorted(hits, key=lambda x: (x[0], x[1], x[2]))
    out = []
    cs, ce, cm = hits[0]
    for s, e, m in hits[1:]:
        if abs(s - cs) <= slop and abs(e - ce) <= slop:
            if (e - s) >= (ce - cs):
                cs, ce, cm = s, e, m
        else:
            out.append((cs, ce, cm))
            cs, ce, cm = s, e, m
    out.append((cs, ce, cm))
    return out

def find_adapter_anywhere(seq, pat, mm_max=1, allow_del=True):
    """
    Scan 'seq' for 'pat' allowing:
      - exact matches (fast path)
      - ≤ mm_max mismatches (Hamming <= mm_max)
      - optionally: exactly 1 deletion in adapter (window len = len(pat)-1), exact
    Returns deduped hits as (start, end, mode) with mode in {"exact","mm","del"}.
    """
    if seq is None:
        return []
    L = len(pat)
    hits = []

    # exact
    start = 0
    while True:
        idx = seq.find(pat, start)
        if idx == -1:
            break
        hits.append((idx, idx + L, "exact"))
        start = idx + 1
    if hits:
        return dedup_hits(hits, slop=0)

    # <= mm_max mismatches
    n = len(seq)
    if n >= L and mm_max >= 0:
        for i in range(n - L + 1):
            if ham(seq[i:i+L], pat) <= mm_max:
                hits.append((i, i + L, "mm"))

    # exactly 1 deletion in adapter
    if allow_del and n >= L - 1:
        w = L - 1
        drop_variants = {pat[:d] + pat[d+1:] for d in range(L)}
        for i in range(n - w + 1):
            if seq[i:i+w] in drop_variants:
                hits.append((i, i + w, "del"))

    return dedup_hits(hits, slop=1)

def scan_for_adapter(seq, adapter, mm_max, allow_del):
    return find_adapter_anywhere(seq, adapter, mm_max=mm_max, allow_del=allow_del)

# ---------------- reporting ----------------

def write_text_log(path, stats, in_bam, out_bam, args):
    total = stats["total"]
    pct = (lambda n: f"{(100.0*n/total):.2f}%" if total else "0.00%")
    with open(path, "w") as f:
        f.write("# Adapter soft-clip filter log\n")
        f.write(f"Input BAM: {in_bam}\n")
        f.write(f"Output BAM: {out_bam}\n")
        f.write(f"Adapter FWD: {args.adapter}\n")
        f.write(f"Adapter RC : {revcomp(args.adapter)}\n")
        f.write(f"Min MAPQ: {args.min_mapq}\n")
        f.write(f"mm_max: {args.mm_max}\n")
        f.write(f"allow_1_deletion: {not args.no_del}\n")
        f.write(f"fullread_classify: {args.fullread_classify}\n")
        f.write(f"extract_n: {args.extract_n}\n")
        f.write(f"extract_out: {args.extract_out}\n\n")
        f.write("Counts:\n")
        for k in (
            "total","kept",
            "removed_unmapped","removed_mapq0_or_below",
            "removed_no_adapter_softclip","removed_adapter_multiple_softclip",
            "removed_no_adapter_anywhere","removed_adapter_multiple_anywhere",
            "removed_adapter_present_not_in_softclip","removed_softclip_multiple_hits",
            "kept_adapter_5p","kept_adapter_3p",
        ):
            f.write(f"  {k}: {stats[k]} ({pct(stats[k])})\n")

def write_tsv(path, stats):
    total = stats["total"]; kept = stats["kept"]
    pT = (lambda n: (100.0*n/total) if total else 0.0)
    pK = (lambda n: (100.0*n/kept)  if kept  else 0.0)
    rows = [
        ("total", stats["total"], pT(stats["total"]), ""),
        ("kept", stats["kept"], pT(stats["kept"]), ""),
        ("removed_unmapped", stats["removed_unmapped"], pT(stats["removed_unmapped"]), ""),
        ("removed_mapq0_or_below", stats["removed_mapq0_or_below"], pT(stats["removed_mapq0_or_below"]), ""),
        ("removed_no_adapter_softclip", stats["removed_no_adapter_softclip"], pT(stats["removed_no_adapter_softclip"]), ""),
        ("removed_adapter_multiple_softclip", stats["removed_adapter_multiple_softclip"], pT(stats["removed_adapter_multiple_softclip"]), ""),
        ("removed_no_adapter_anywhere", stats["removed_no_adapter_anywhere"], pT(stats["removed_no_adapter_anywhere"]), ""),
        ("removed_adapter_multiple_anywhere", stats["removed_adapter_multiple_anywhere"], pT(stats["removed_adapter_multiple_anywhere"]), ""),
        ("removed_adapter_present_not_in_softclip", stats["removed_adapter_present_not_in_softclip"], pT(stats["removed_adapter_present_not_in_softclip"]), ""),
        ("removed_softclip_multiple_hits", stats["removed_softclip_multiple_hits"], pT(stats["removed_softclip_multiple_hits"]), ""),
        ("kept_adapter_5p", stats["kept_adapter_5p"], pT(stats["kept_adapter_5p"]), pK(stats["kept_adapter_5p"])),
        ("kept_adapter_3p", stats["kept_adapter_3p"], pT(stats["kept_adapter_3p"]), pK(stats["kept_adapter_3p"])),
    ]
    with open(path, "w") as f:
        f.write("metric\tcount\tpct_of_total\tpct_of_kept\n")
        for metric, count, p1, p2 in rows:
            if p2 == "":
                f.write(f"{metric}\t{count}\t{p1:.4f}\t\n")
            else:
                f.write(f"{metric}\t{count}\t{p1:.4f}\t{p2:.4f}\n")

def write_extract_header(tsv_fp):
    if tsv_fp:
        tsv_fp.write("read_id\tside\thit_mode\tadapter_start\tadapter_end\textracted\n")

def emit_extraction(tsv_fp, read_id, side, mode, s, e, seq):
    if tsv_fp:
        tsv_fp.write(f"{read_id}\t{side}\t{mode}\t{s}\t{e}\t{seq}\n")

# ---------------- core ----------------

def process_bam(in_bam, out_bam, adapter_fwd, min_mapq,
                mm_max, allow_del, fullread_classify,
                extract_n, extract_out):
    stats = Counter()
    adapter_rc = revcomp(adapter_fwd)

    tsv_fp = open(extract_out, "w") if (extract_out and extract_n > 0) else None
    if tsv_fp:
        write_extract_header(tsv_fp)

    with pysam.AlignmentFile(in_bam, "rb") as inf, \
         pysam.AlignmentFile(out_bam, "wb", header=inf.header) as outf:

        for r in inf:
            stats["total"] += 1

            if r.is_unmapped:
                stats["removed_unmapped"] += 1
                continue

            if r.mapping_quality is None or r.mapping_quality <= 0 or r.mapping_quality < min_mapq:
                stats["removed_mapq0_or_below"] += 1
                continue

            sc5, sc3 = get_softclips(r)

            # Search soft-clips (anywhere within the clip)
            hits5 = scan_for_adapter(sc5, adapter_rc, mm_max, allow_del)
            hits3 = scan_for_adapter(sc3, adapter_fwd, mm_max, allow_del)
            softclip_total = len(hits5) + len(hits3)

            # Keep only if exactly one hit across soft-clips
            if softclip_total == 1:
                if len(hits5) == 1:
                    stats["kept_adapter_5p"] += 1
                    outf.write(r)
                    stats["kept"] += 1

                    # 5' (RC) extraction: [s : e+N] within sc5, then revcomp
                    if extract_n > 0 and sc5:
                        s, e, mode = hits5[0]
                        end = min(len(sc5), e + extract_n)
                        tmp = sc5[s:end]
                        extracted = revcomp(tmp)
                        emit_extraction(tsv_fp, r.query_name, "5p", mode, s, e, extracted)

                else:
                    stats["kept_adapter_3p"] += 1
                    outf.write(r)
                    stats["kept"] += 1

                    # 3' (FWD) extraction: [s-N : e] within sc3, as-is
                    if extract_n > 0 and sc3:
                        s, e, mode = hits3[0]
                        start = max(0, s - extract_n)
                        extracted = sc3[start:e]
                        emit_extraction(tsv_fp, r.query_name, "3p", mode, s, e, extracted)

                continue

            # Rejected → classify
            if not fullread_classify:
                if softclip_total == 0:
                    stats["removed_no_adapter_softclip"] += 1
                else:
                    stats["removed_adapter_multiple_softclip"] += 1
                continue

            # Full-read classification (both orientations)
            qseq = r.query_sequence or ""
            hits_any = []
            if qseq:
                hits_any += scan_for_adapter(qseq, adapter_fwd, mm_max, allow_del)
                hits_any += scan_for_adapter(qseq, adapter_rc,  mm_max, allow_del)
                hits_any = dedup_hits(hits_any, slop=1)
            anywhere_total = len(hits_any)

            if anywhere_total == 0:
                stats["removed_no_adapter_anywhere"] += 1
                continue
            if anywhere_total >= 2:
                stats["removed_adapter_multiple_anywhere"] += 1
                continue

            # anywhere_total == 1 but softclip_total != 1
            if softclip_total == 0:
                stats["removed_adapter_present_not_in_softclip"] += 1
            else:
                stats["removed_softclip_multiple_hits"] += 1
            continue

    if tsv_fp:
        tsv_fp.close()
    return stats

def main():
    ap = argparse.ArgumentParser(
        description="Filter BAM by adapter-in-softclips; optional full-read classification; extract adjacent window into TSV."
    )
    ap.add_argument("-i", "--input-bam", required=True)
    ap.add_argument("-o", "--output-bam", required=True)
    ap.add_argument("--adapter", default=ADAPTER_DEFAULT, help="3' adapter (forward orientation).")
    ap.add_argument("--min-mapq", type=int, default=1, help="Minimum MAPQ to keep (default 1 removes MAPQ=0).")
    ap.add_argument("--mm-max", type=int, default=1, help="Maximum mismatches allowed (default 1; use 2 to relax).")
    ap.add_argument("--no-del", action="store_true", help="Disable the 1-deletion allowance.")
    ap.add_argument("--fullread-classify", action="store_true",
                    help="Also scan the full read (both orientations) to classify rejects.")
    # Extraction (TSV)
    ap.add_argument("--extract-n", type=int, default=0, help="Window size N for adjacent extraction (0 disables).")
    ap.add_argument("--extract-out", default=None, help="Write extracted sequences to this TSV (requires --extract-n > 0).")
    ap.add_argument("--log", default=None)
    ap.add_argument("--tsv", default=None)
    args = ap.parse_args()

    log_path = args.log if args.log else args.output_bam + ".log.txt"
    tsv_path = args.tsv if args.tsv else args.output_bam + ".stats.tsv"

    stats = process_bam(
        args.input_bam, args.output_bam, args.adapter, args.min_mapq,
        mm_max=args.mm_max, allow_del=(not args.no_del),
        fullread_classify=args.fullread_classify,
        extract_n=args.extract_n, extract_out=args.extract_out
    )
    write_text_log(log_path, stats, args.input_bam, args.output_bam, args)
    write_tsv(tsv_path, stats)

if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        try:
            sys.stdout.close()
        except Exception:
            pass
        try:
            sys.stderr.close()
        except Exception:
            pass
