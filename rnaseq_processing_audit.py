#!/usr/bin/env python3
"""
rnaseq_processing_audit.py  --  post-hoc "were all samples processed?" audit

Reads ONLY the outputs an already-finished run left on disk (no re-run, no FASTQ).
It reconstructs the audit trail every sample must leave behind --
  sample_sheet -> trimmed pair -> BAM(+index) -> HISAT2 log block -> counts column -> summary column
-- and flags any sample that is missing at any stage, plus the sheet-hygiene and
column-labelling problems (S1 / R1 / S2 / S3) raised in the code review.

Stdlib only. Usage:
    python3 rnaseq_processing_audit.py [PROJECT_DIR]
    python3 rnaseq_processing_audit.py [PROJECT_DIR] --csv audit_matrix.csv

PROJECT_DIR is the base under which trimmed/ mapped/ pipeline.log qc/ and (by default)
sample_sheet.tsv + counts/counts.txt are looked for. Any of those can be pointed
elsewhere with the override flags below -- useful when the counts matrix was regenerated
in a separate folder, or the run is spread across directories:
    --sheet PATH     --mapped DIR     --trimmed DIR
    --counts FILE    --summary FILE   --log FILE     --qc DIR
If no sample sheet is available the expected sample set is derived from mapped/*.bam
instead (this is the right reference when counts were regenerated from existing BAMs).
Any stage whose folder/file is absent is SKIPPED (not failed).

Exit code 0 = every expected sample traced through every AVAILABLE stage (a residual S3
advisory about MultiQC coverage is informational and does NOT fail the run); 1 = at least
one CRITICAL finding or an incomplete stage; 2 = neither a sample sheet nor a mapped/ dir
could be found (nothing to audit against). Usable in CI.
"""

import os, re, sys, csv, argparse

# ---- expected pipeline layout (mirrors run_rnaseq_pipeline_010726.sh) -------
def paths(proj, args):
    """Resolve every input path: an explicit --flag wins, else fall back to a
    sensible default under PROJECT_DIR. counts default tries counts/counts.txt then
    counts.txt; summary defaults to <counts>.summary."""
    def d(*p): return os.path.join(proj, *p)
    counts = args.counts
    if counts is None:
        counts = d("counts", "counts.txt")
        if not os.path.isfile(counts) and os.path.isfile(d("counts.txt")):
            counts = d("counts.txt")
    summary = args.summary if args.summary else (counts + ".summary")
    return dict(
        sheet   = args.sheet   if args.sheet   else d("sample_sheet.tsv"),
        trimmed = args.trimmed if args.trimmed else d("trimmed"),
        mapped  = args.mapped  if args.mapped  else d("mapped"),
        counts  = counts,
        summary = summary,
        log     = args.log     if args.log     else d("pipeline.log"),
        qc      = args.qc      if args.qc      else d("qc"),
    )

C_RED, C_AMBER, C_GREEN, C_DIM, C_OFF = "\033[91m", "\033[93m", "\033[92m", "\033[2m", "\033[0m"
def _c(s, c):  # colour only when writing to a terminal
    return f"{c}{s}{C_OFF}" if sys.stdout.isatty() else s
findings = []  # (severity, code, message)
def flag(sev, code, msg): findings.append((sev, code, msg))

# ---- sample sheet ----------------------------------------------------------
def read_sheet(path):
    """Return (samples[list], hygiene[dict]). Robust to CRLF / missing final newline."""
    with open(path, "rb") as fh:
        raw = fh.read()
    hy = dict(
        crlf            = b"\r\n" in raw,
        lone_cr         = bool(re.search(rb"\r(?!\n)", raw)),
        trailing_nl     = raw.endswith(b"\n") or raw.endswith(b"\r"),
        bytes           = len(raw),
    )
    text = raw.decode("utf-8-sig", errors="replace").replace("\r\n", "\n").replace("\r", "\n")
    lines = [ln for ln in text.split("\n") if ln.strip() != ""]
    hdr = lines[0].split("\t")
    rows = [ln.split("\t") for ln in lines[1:]]
    samples, badcols = [], []
    for i, r in enumerate(rows, start=2):
        if len(r) != len(hdr):
            badcols.append((i, len(r)))
        samples.append(r[0].strip() if r else "")
    hy["n_header_cols"]  = len(hdr)
    hy["n_data_rows"]    = len(rows)
    hy["ragged_rows"]    = badcols
    hy["dupes"]          = sorted({s for s in samples if samples.count(s) > 1})
    hy["header"]         = hdr
    return samples, hy

# ---- log parsing -----------------------------------------------------------
def parse_log(path):
    """Per-sample HISAT2 overall alignment rate + the 'Finished' roll-call + completion marker."""
    started, finished, rates = [], [], {}
    completed = False
    if not os.path.isfile(path):
        return started, finished, rates, completed, False
    cur = None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            m = re.search(r"Processing sample \d+ of \d+:\s*(\S+)", ln)
            if m:
                cur = m.group(1); started.append(cur); continue
            m = re.search(r"([\d.]+)%\s+overall alignment rate", ln)
            if m and cur is not None:
                rates[cur] = float(m.group(1)); continue
            m = re.search(r"Finished processing:\s*(\S+)", ln)
            if m:
                finished.append(m.group(1)); continue
            if "PIPELINE COMPLETED SUCCESSFULLY" in ln:
                completed = True
    return started, finished, rates, completed, True

# ---- featureCounts header --------------------------------------------------
def bam_cols_to_samples(header_fields, lead_cols):
    """Turn '.../mapped/S01_x.bam' column headers into sample names."""
    out = []
    for f in header_fields[lead_cols:]:
        base = os.path.basename(f.strip())
        base = re.sub(r"\.bam$", "", base, flags=re.IGNORECASE)
        out.append(base)
    return out

def read_counts_header(path):
    if not os.path.isfile(path): return None, None
    with open(path, encoding="utf-8", errors="replace") as fh:
        first = fh.readline()
        if first.startswith("#"):           # featureCounts comment line
            header = fh.readline().rstrip("\n").split("\t")
            had_comment = True
        else:
            header = first.rstrip("\n").split("\t"); had_comment = False
    # standard featureCounts lead columns: Geneid Chr Start End Strand Length
    lead = 6 if header[:6] == ["Geneid","Chr","Start","End","Strand","Length"] else \
           next((i for i,h in enumerate(header) if h.endswith(".bam")), len(header))
    return bam_cols_to_samples(header, lead), had_comment

def read_summary(path):
    """counts.txt.summary -> {sample: {status: int}} keyed by bam-derived sample name."""
    if not os.path.isfile(path): return {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
    cols = bam_cols_to_samples(rows[0], 1)
    data = {c: {} for c in cols}
    for r in rows[1:]:
        status = r[0]
        for c, v in zip(cols, r[1:]):
            try: data[c][status] = int(v)
            except ValueError: pass
    return data

# ---- main ------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("project_dir", nargs="?", default="/mnt/d/OxoMech")
    ap.add_argument("--csv", help="write per-sample audit matrix to this CSV")
    ap.add_argument("--min-align", type=float, default=70.0, help="warn below this overall alignment %% (default 70)")
    ap.add_argument("--min-assign", type=float, default=50.0, help="warn below this featureCounts assignment %% (default 50)")
    # per-stage path overrides (for runs split across folders / counts regenerated elsewhere)
    ap.add_argument("--sheet",   help="path to sample_sheet.tsv (default PROJECT_DIR/sample_sheet.tsv)")
    ap.add_argument("--mapped",  help="mapped/ dir with *.bam (default PROJECT_DIR/mapped)")
    ap.add_argument("--trimmed", help="trimmed/ dir with *_R1/R2.fastq.gz (default PROJECT_DIR/trimmed)")
    ap.add_argument("--counts",  help="featureCounts counts.txt (default PROJECT_DIR/counts/counts.txt then PROJECT_DIR/counts.txt)")
    ap.add_argument("--summary", help="counts.txt.summary (default <counts>.summary)")
    ap.add_argument("--log",     help="pipeline.log (default PROJECT_DIR/pipeline.log)")
    ap.add_argument("--qc",      help="qc/ dir with MultiQC html (default PROJECT_DIR/qc)")
    args = ap.parse_args()
    P = paths(args.project_dir, args)

    print(f"\n{'='*74}\n  RNA-seq processing audit  --  {args.project_dir}\n{'='*74}")

    # BAM-derived sample list (used either to cross-check the sheet, or AS the
    # reference sample set when no sheet is available).
    def bams_in(d):
        if not os.path.isdir(d): return []
        return sorted(re.sub(r"\.bam$", "", f, flags=re.I)
                      for f in os.listdir(d) if f.lower().endswith(".bam"))
    bam_samples = bams_in(P["mapped"])

    have_sheet = os.path.isfile(P["sheet"])
    if have_sheet:
        samples, hy = read_sheet(P["sheet"])
        ref_src = "sample sheet"
    elif bam_samples:
        samples, hy = bam_samples, None
        ref_src = "mapped/*.bam"
    else:
        print(_c(f"FATAL: no sample sheet ({P['sheet']}) AND no BAMs in {P['mapped']} -- "
                 f"nothing to audit against. Point --sheet and/or --mapped at the right folders.", C_RED))
        sys.exit(2)
    n = len(samples)

    # ---------- Section A: sample-sheet hygiene (S1) ----------
    print(f"\n[A] SAMPLE SHEET HYGIENE")
    if not have_sheet:
        flag("WARN","S1", f"No sample_sheet.tsv at {P['sheet']} -- sheet hygiene (CRLF / trailing-newline) "
                          f"could NOT be checked. Expected sample set derived from {len(bam_samples)} BAM(s) "
                          f"in {P['mapped']} instead. Pass --sheet to enable S1 checks.")
        print(_c(f"    no sheet found at {P['sheet']}", C_AMBER))
        print(f"    using {len(bam_samples)} BAM name(s) as the reference sample set instead")
    else:
        print(f"    path           : {P['sheet']}")
        print(f"    header columns : {hy['n_header_cols']}  ->  {hy['header']}")
        print(f"    data rows      : {hy['n_data_rows']}")
        if hy["crlf"]:
            flag("CRIT","S1", "Sheet has CRLF (Windows) line endings -- trailing '\\r' corrupts the last field "
                              "(R2 path) of every row. Fix: dos2unix sample_sheet.tsv")
            print(_c("    line endings   : CRLF  <-- CRITICAL (run dos2unix)", C_RED))
        else:
            print(_c("    line endings   : LF (unix)  OK", C_GREEN))
        if not hy["trailing_nl"]:
            flag("CRIT","S1", "Sheet has NO trailing newline -- bash `while read` drops the LAST sample silently. "
                              "Fix: printf '\\n' >> sample_sheet.tsv")
            print(_c("    trailing NL    : MISSING  <-- CRITICAL (last sample skipped by `while read`)", C_RED))
        else:
            print(_c("    trailing NL    : present  OK", C_GREEN))
        if hy["ragged_rows"]:
            flag("CRIT","S1", f"{len(hy['ragged_rows'])} row(s) have wrong column count: {hy['ragged_rows']}")
            print(_c(f"    ragged rows    : {hy['ragged_rows']}  <-- CRITICAL", C_RED))
        if hy["dupes"]:
            flag("CRIT","S1", f"Duplicate sample names in sheet: {hy['dupes']}")
            print(_c(f"    duplicate IDs  : {hy['dupes']}  <-- CRITICAL", C_RED))

    # ---------- gather per-stage evidence (skip stages whose folder/file is absent) ----------
    trimmed_r1 = {s for s in samples if os.path.isfile(os.path.join(P["trimmed"], f"{s}_R1.fastq.gz"))} if os.path.isdir(P["trimmed"]) else set()
    trimmed_r2 = {s for s in samples if os.path.isfile(os.path.join(P["trimmed"], f"{s}_R2.fastq.gz"))} if os.path.isdir(P["trimmed"]) else set()
    bam        = {s for s in samples if os.path.isfile(os.path.join(P["mapped"],  f"{s}.bam"))} if os.path.isdir(P["mapped"]) else set()
    bai        = {s for s in samples if os.path.isfile(os.path.join(P["mapped"],  f"{s}.bam.bai"))} if os.path.isdir(P["mapped"]) else set()
    # A stage is only INFORMATIVE if at least one matching file was found. Zero matches
    # almost always means the intermediate was cleaned up (trimmed FASTQs routinely are)
    # or named differently -- NOT that every sample was dropped. Treat empty as unpopulated.
    have_trimmed = len(trimmed_r1) > 0 or len(trimmed_r2) > 0
    have_mapped  = len(bam) > 0
    started, finished, rates, completed, have_log = parse_log(P["log"])
    finished_set = set(finished)
    counts_cols, had_comment = read_counts_header(P["counts"])
    counts_set  = set(counts_cols) if counts_cols else set()
    summ        = read_summary(P["summary"])
    summ_set    = set(summ)
    print(f"\n    resolved inputs (reference sample set = {ref_src}, n={n}):")
    for k in ["trimmed","mapped","log","counts","summary","qc"]:
        ok = (os.path.isdir(P[k]) if k in ("trimmed","mapped","qc") else os.path.isfile(P[k]))
        print(f"      {k:8s}: {P[k]}  {'' if ok else _c('[absent -- skipped]', C_AMBER)}")

    # ---------- Section B: per-sample cross-stage matrix ----------
    stages = ["trimR1","trimR2","BAM","BAI","log_done","align%","counts","assign%"]
    print(f"\n[B] PER-SAMPLE AUDIT MATRIX   (v=present  .=MISSING)")
    print("    " + "sample".ljust(26) + "  " + "  ".join(h.ljust(8) for h in stages))
    matrix = []
    for s in samples:
        arate = rates.get(s)
        assigned = summ.get(s, {}).get("Assigned")
        tot = sum(summ.get(s, {}).values()) if s in summ else None
        apct = (100.0*assigned/tot) if (assigned is not None and tot) else None
        row = dict(
            sample=s,
            trimR1 = s in trimmed_r1, trimR2 = s in trimmed_r2,
            BAM = s in bam, BAI = s in bai,
            log_done = s in finished_set,
            align = arate, counts = s in counts_set, assign = apct,
        )
        matrix.append(row)
        def mark(b, avail=True):
            if not avail: return "-".ljust(8)                       # stage not present -> not evaluated
            return _c("v".ljust(8), C_GREEN) if b else _c(".".ljust(8), C_RED)
        def rate_cell(val, avail):
            if not avail: return "-".ljust(8)
            return f"{val:.1f}".ljust(8) if val is not None else _c(".".ljust(8), C_RED)
        cells = [mark(row["trimR1"], have_trimmed), mark(row["trimR2"], have_trimmed),
                 mark(row["BAM"], have_mapped), mark(row["BAI"], have_mapped),
                 mark(row["log_done"], have_log),
                 rate_cell(arate, have_log),
                 mark(row["counts"], counts_cols is not None),
                 rate_cell(apct, bool(summ))]
        print("    " + s.ljust(26) + "  " + "  ".join(cells))

    # ---------- Section C: completeness verdict (S1 core) ----------
    def missing(stage_set): return [s for s in samples if s not in stage_set]
    print(f"\n[C] COMPLETENESS  (reference sample set = {ref_src}, n={n})")
    # (label, evidence set, is-this-stage-available). Absent stages are SKIPPED, not failed.
    per_stage = [("trimmed R1", trimmed_r1, have_trimmed), ("trimmed R2", trimmed_r2, have_trimmed),
                 ("BAM", bam, have_mapped), ("BAM index", bai, have_mapped),
                 ("HISAT2 'Finished' log", finished_set, have_log),
                 ("counts.txt column", counts_set, counts_cols is not None),
                 ("summary column", summ_set, bool(summ))]
    all_ok = True
    for name, st, avail in per_stage:
        if not avail:
            print(_c(f"    {name.ljust(24)}: -- stage not present, skipped", C_AMBER)); continue
        miss = missing(st)
        if miss:
            all_ok = False
            flag("CRIT","S1", f"{len(miss)} sample(s) missing at stage '{name}': {miss}")
            print(_c(f"    {name.ljust(24)}: {len(st & set(samples))}/{n}  MISSING {miss}", C_RED))
        else:
            print(_c(f"    {name.ljust(24)}: {n}/{n}  OK", C_GREEN))

    # ---------- Section D: counts columns vs reference set (R1) ----------
    ref = ref_src  # "sample sheet" or "mapped/*.bam"
    print(f"\n[D] COUNTS-vs-{ref.upper()} LABEL CHECK (R1: intersect() silently drops mismatches)")
    if counts_cols is None:
        flag("CRIT","S1", f"counts.txt not found at {P['counts']} -- cannot verify columns. "
                          f"Pass --counts if it lives elsewhere.")
        print(_c(f"    counts.txt not found at {P['counts']}  <-- CRITICAL", C_RED))
    else:
        in_counts_not_ref = [c for c in counts_cols if c not in samples]
        in_ref_not_counts = [s for s in samples if s not in counts_set]
        print(f"    counts.txt columns : {len(counts_cols)}   (reference {ref}: {n})")
        print(f"    featureCounts '#' comment line present (R5: skip=1 assumption): "
              + (_c("yes  OK", C_GREEN) if had_comment else _c("NO  <-- skip=1 would drop the header!", C_AMBER)))
        if in_ref_not_counts:
            flag("CRIT","R1", f"Sample(s) in {ref} but NOT a counts.txt column (would be silently dropped by "
                              f"intersect() in the DESeq2 script): {in_ref_not_counts}")
            print(_c(f"    in {ref}, not in counts : {in_ref_not_counts}  <-- CRITICAL (R1)", C_RED))
        if in_counts_not_ref:
            sev = "WARN" if have_sheet else "CRIT"
            flag(sev,"R1", f"counts.txt column(s) with no matching {ref} entry: {in_counts_not_ref}")
            print(_c(f"    in counts, not in {ref} : {in_counts_not_ref}  <-- {sev}", C_RED if sev=="CRIT" else C_AMBER))
        if not in_ref_not_counts and not in_counts_not_ref:
            order_ok = counts_cols == [s for s in samples if s in counts_set]
            print(_c(f"    set match               : identical ({len(counts_cols)}/{n})  OK", C_GREEN))
            if have_sheet:
                print(_c(f"    column order == sheet    : {'yes  OK' if order_ok else 'NO (reorder in R by name)'}",
                         C_GREEN if order_ok else C_AMBER))
                if not order_ok:
                    flag("WARN","R1","counts.txt column order differs from sheet order -- safe only if you key "
                                     "colData by name (DESeqDataSetFromMatrix guards this), not by position.")

    # ---------- Section E: alignment & assignment rates (S2 / S3) ----------
    print(f"\n[E] ALIGNMENT & ASSIGNMENT RATES")
    if not have_log:
        flag("WARN","S3", "pipeline.log not found -- cannot recover per-sample HISAT2 rates "
                          "(they are not written as standalone files; MultiQC never sees them either).")
        print(_c("    pipeline.log not found  <-- WARN (S3)", C_AMBER))
    else:
        low_align = [s for s in samples if rates.get(s) is not None and rates[s] < args.min_align]
        got = [rates[s] for s in samples if s in rates]
        if got:
            print(f"    overall alignment %  : min {min(got):.1f}  median {sorted(got)[len(got)//2]:.1f}  max {max(got):.1f}")
        if low_align:
            flag("WARN","QC", f"Low alignment (<{args.min_align:.0f}%): "
                              + ", ".join(f'{s}={rates[s]:.1f}%' for s in low_align))
            print(_c(f"    low-alignment samples: {low_align}  <-- WARN", C_AMBER))
    if summ:
        apcts = {}
        for s in samples:
            if s in summ and sum(summ[s].values()):
                apcts[s] = 100.0*summ[s].get("Assigned",0)/sum(summ[s].values())
        if apcts:
            vals = list(apcts.values())
            print(f"    featureCounts assign%: min {min(vals):.1f}  median {sorted(vals)[len(vals)//2]:.1f}  max {max(vals):.1f}")
            low_assign = [s for s in apcts if apcts[s] < args.min_assign]
            if low_assign:
                flag("WARN","S2", f"Low assignment (<{args.min_assign:.0f}%) with -s 2 (reverse-stranded). "
                                  f"If unexpectedly low the library may not be reverse-stranded -- confirm strandedness "
                                  f"with RSeQC infer_experiment.py: {[f'{s}={apcts[s]:.1f}%' for s in low_assign]}")
                print(_c(f"    low-assignment (S2)  : {low_assign}  <-- WARN (verify strandedness)", C_AMBER))
    else:
        print(_c("    counts.txt.summary not found -- assignment rates unavailable", C_AMBER))

    # ---------- Section F: QC coverage (S3) & completion marker ----------
    print(f"\n[F] QC OUTPUTS & COMPLETION")
    mqc = [f for f in (os.listdir(P["qc"]) if os.path.isdir(P["qc"]) else []) if f.endswith(".html")]
    print(f"    MultiQC reports in qc/ : {mqc if mqc else 'none'}")
    flag("WARN","S3", "MultiQC was pointed only at FastQC dirs; HISAT2 rates (in pipeline.log) and the "
                      "featureCounts .summary are NOT in the report -- this script recovers them, but consider "
                      "adding mapped/ logs + counts.txt.summary to a multiqc run for a single QC dashboard.")
    if have_log:
        if completed:
            print(_c("    'PIPELINE COMPLETED SUCCESSFULLY' marker: present  OK", C_GREEN))
        else:
            flag("WARN","S1","pipeline.log has NO 'PIPELINE COMPLETED SUCCESSFULLY' line -- the run may have "
                             "died mid-loop (set -e aborts on first error).")
            print(_c("    'PIPELINE COMPLETED SUCCESSFULLY' marker: ABSENT  <-- WARN (run may have aborted)", C_AMBER))

    # ---------- CSV ----------
    if args.csv:
        with open(args.csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["sample","trimR1","trimR2","BAM","BAI","log_finished","align_pct","in_counts","assign_pct"])
            for r in matrix:
                w.writerow([r["sample"], int(r["trimR1"]), int(r["trimR2"]), int(r["BAM"]), int(r["BAI"]),
                            int(r["log_done"]),
                            f"{r['align']:.2f}" if r["align"] is not None else "",
                            int(r["counts"]),
                            f"{r['assign']:.2f}" if r["assign"] is not None else ""])
        print(f"\n    per-sample matrix written -> {args.csv}")

    # ---------- verdict ----------
    crits = [f for f in findings if f[0]=="CRIT"]
    warns = [f for f in findings if f[0]=="WARN"]
    print(f"\n{'='*74}\n  VERDICT")
    if not crits and all_ok:
        print(_c(f"  PASS -- all {n} samples traced through every stage; no critical issues.", C_GREEN))
    else:
        print(_c(f"  FAIL -- {len(crits)} CRITICAL, {len(warns)} review item(s).", C_RED))
    for sev, code, msg in sorted(findings, key=lambda x: {"CRIT":0,"WARN":1}[x[0]]):
        tag = _c(f"[{sev} {code}]", C_RED if sev=="CRIT" else C_AMBER)
        print(f"  {tag} {msg}")
    print("="*74 + "\n")
    sys.exit(1 if (crits or not all_ok) else 0)

if __name__ == "__main__":
    main()
