#!/usr/bin/env python
import os
import subprocess
import sys
import warnings
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import Counter
from datetime import datetime, timezone 

from Bio.Seq import Seq
from pyfaidx import Fasta

import tempfile

# Pick a base tmp dir that is *meant* for the current job/user
TMP_BASE = os.environ.get("TMPDIR", "/tmp")

# Create a unique, private temp dir for THIS run of the script
LOCAL_CWD = tempfile.mkdtemp(prefix="chewbbaca_runs_", dir=TMP_BASE)

# (optional) make sure permissions are sane even with umask weirdness
os.chmod(LOCAL_CWD, 0o700)

# Define valid bases and codon sets (uppercase)
VALID_BASES = set(("A", "T", "C", "G"))
START_CODONS = {"ATG", "TTG", "GTG", "CTG", "ATA", "ATT"}
STOP_CODONS  = {"TAA", "TAG", "TGA"}


def is_valid_cds(seq: str) -> bool:
    return (
        seq[:3] in START_CODONS and
        seq[-3:] in STOP_CODONS and
        len(seq) % 3 == 0
    )


def has_internal_stop(seq: str) -> bool:
    for i in range(3, len(seq) - 3, 3):
        if seq[i:i+3] in STOP_CODONS:
            return True
    return False


def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def orient_and_frame_fix(
    seq: str,
    qstart: int,
    qend: int,
    genome: Fasta,
    seq_id: str,
    header_raw: str,
    strand: int,
    max_stop_extend: int
):
    """
    Strand-aware orientation + frame correction + optional in-frame stop extension.
    Returns (fixed_seq, new_start, new_end, stats).
    If it cannot be fixed, returns ("", qstart, qend, stats) to signal caller to skip.
    """
    stats = Counter()

    # 1) Orientation from BLAST strand
    if strand == -1:
        seq = reverse_complement(seq)
        stats["reoriented"] += 1

    seq = seq.upper()

    # 2) Frame correction by padding (no trimming fallback; skip if cannot restore)
    length = len(seq)
    mod = length % 3
    if mod != 0:
        needed = 3 - mod
        adjusted = False
        contig_len = len(genome[seq_id])

        # a) Extend 3' end by needed bases
        if strand == 1:
            if qend + needed <= contig_len:
                cand = genome[seq_id][qstart:qend + needed].seq.upper()
                if len(cand) % 3 == 0:
                    seq = cand
                    qend += needed
                    stats["padded_stop"] += 1
                    adjusted = True
        else:
            if qstart - needed >= 0:
                cand = genome[seq_id][qstart - needed:qend].seq.upper()
                if len(cand) % 3 == 0:
                    seq = reverse_complement(cand)
                    qstart -= needed
                    stats["padded_stop"] += 1
                    adjusted = True

        # b) Extend 5' end by needed bases
        if not adjusted:
            if strand == 1 and qstart - needed >= 0:
                cand = genome[seq_id][qstart - needed:qend].seq.upper()
                if len(cand) % 3 == 0:
                    seq = cand
                    qstart -= needed
                    stats["padded_start"] += 1
                    adjusted = True
            elif strand == -1 and qend + needed <= contig_len:
                cand = genome[seq_id][qstart:qend + needed].seq.upper()
                if len(cand) % 3 == 0:
                    seq = reverse_complement(cand)
                    qend += needed
                    stats["padded_start"] += 1
                    adjusted = True

    if len(seq) % 3 != 0:
        warnings.warn(f"{header_raw}: out of frame after padding; skipping.", RuntimeWarning)
        print(f"[WARNING] {header_raw}: out of frame after padding; skipping.", file=sys.stderr)
        return "", qstart, qend, stats

    # 3) If no terminal stop, extend in-frame up to max_stop_extend (bp)
    if seq[-3:] not in STOP_CODONS:
        contig_len = len(genome[seq_id])
        extended = False

        for i in range(3, max_stop_extend + 1, 3):
            if strand == 1:
                if qend + i > contig_len:
                    break
                cand = genome[seq_id][qstart:qend + i].seq.upper()
            else:
                if qstart - i < 0:
                    break
                raw = genome[seq_id][qstart - i:qend].seq.upper()
                cand = reverse_complement(raw)

            if cand[-3:] in STOP_CODONS and len(cand) % 3 == 0:
                seq = cand
                if strand == 1:
                    qend += i
                else:
                    qstart -= i
                stats["stop_extended"] = i
                extended = True
                break

        if not extended:
            msg = f"{header_raw}: no terminal stop within {max_stop_extend} bp 3' window; skipping."
            warnings.warn(msg, RuntimeWarning)
            print(f"[WARNING] {msg}", file=sys.stderr)
            return "", qstart, qend, stats

    return seq, qstart, qend, stats


def extract_subsequences(genome: Fasta, alleles: list, max_stop_extend: int) -> list:
    extracted = []

    for allele in alleles:
        seq_id, qstart, qend, saccver, strand = allele
        header_raw = f">{seq_id}_{saccver}_{qstart+1}_{qend}"

        subseq = genome[seq_id][qstart:qend].seq

        if len(subseq) != (qend - qstart):
            msg = f"{header_raw}: slice length {len(subseq)} != qend-qstart ({qend - qstart}); skipping."
            warnings.warn(msg, RuntimeWarning)
            print(f"[WARNING] {msg}", file=sys.stderr)
            continue

        invalid_positions = [(i, b) for i, b in enumerate(subseq) if b.upper() not in VALID_BASES]
        if invalid_positions:
            for pos, base in invalid_positions[:20]:
                msg = f"{header_raw}: invalid base '{base}' at position {pos+1}"
                warnings.warn(msg, RuntimeWarning)
                print(f"[WARNING] {msg}", file=sys.stderr)
            if len(invalid_positions) > 20:
                msg = f"{header_raw}: ... and {len(invalid_positions)-20} more invalid base(s)"
                warnings.warn(msg, RuntimeWarning)
                print(f"[WARNING] {msg}", file=sys.stderr)
            continue

        fixed, new_start, new_end, _stats = orient_and_frame_fix(
            subseq, qstart, qend, genome, seq_id, header_raw, strand, max_stop_extend
        )
        if not fixed:
            continue

        if has_internal_stop(fixed):
            msg = f"{header_raw}: contains internal in-frame stop codon; skipping."
            warnings.warn(msg, RuntimeWarning)
            print(f"[WARNING] {msg}", file=sys.stderr)
            continue

        if is_valid_cds(fixed):
            final_hdr = f">{seq_id}_{saccver}_{new_start+1}_{new_end}"
            extracted.append((final_hdr, fixed))
        else:
            msg = f"{header_raw}: invalid CDS; skipping."
            warnings.warn(msg, RuntimeWarning)
            print(f"[WARNING] {msg}", file=sys.stderr)

    return extracted


def parse_blast_output(blast_output: str, genome: Fasta, min_cov_ratio: float) -> list:
    perfect_hits = []
    fallback_strict = None
    fallback_strict_key = None
    fallback_any = None
    fallback_any_key = None

    with open(blast_output) as fh:
        for line in fh:
            cols = line.strip().split("\t")
            if len(cols) < 13:
                warnings.warn(f"Malformed BLAST line (expected 13 cols): {line.strip()}", RuntimeWarning)
                continue

            rec = {
                "qaccver": cols[0],
                "saccver": cols[1],
                "slen": int(cols[2]),
                "pident": float(cols[3]),
                "length": int(cols[4]),
                "mismatch": int(cols[5]),
                "gapopen": int(cols[6]),
                "qstart": int(cols[7]),
                "qend": int(cols[8]),
                "sstart": int(cols[9]),
                "send": int(cols[10]),
                "evalue": float(cols[11]),
                "bitscore": float(cols[12]),
            }

            qstart, qend = rec["qstart"], rec["qend"]
            query_span = abs(qend - qstart) + 1
            no_gaps = (rec["gapopen"] == 0)

            if rec["pident"] == 100.0 and no_gaps and query_span == rec["slen"]:
                perfect_hits.append(rec)
            else:
                rank_key = (no_gaps, query_span, rec["pident"], rec["bitscore"])

                if fallback_any_key is None or rank_key > fallback_any_key:
                    fallback_any = rec
                    fallback_any_key = rank_key

                if query_span >= min_cov_ratio * rec["slen"]:
                    if fallback_strict_key is None or rank_key > fallback_strict_key:
                        fallback_strict = rec
                        fallback_strict_key = rank_key

    # de-duplicate perfect hits that share qstart or qend on same contig + locus
    if len(perfect_hits) > 1:
        def span(r): return abs(r["qend"] - r["qstart"]) + 1
        perfect_hits.sort(key=lambda r: (span(r), r["bitscore"]), reverse=True)

        kept = []
        seen_start = set()
        seen_end = set()
        for r in perfect_hits:
            kS = (r["qaccver"], r["saccver"], r["qstart"])
            kE = (r["qaccver"], r["saccver"], r["qend"])
            if kS in seen_start or kE in seen_end:
                continue
            kept.append(r)
            seen_start.add(kS)
            seen_end.add(kE)
        perfect_hits = kept

    hits = perfect_hits if perfect_hits else (
        [fallback_strict] if fallback_strict else ([fallback_any] if fallback_any else [])
    )

    alleles = []
    for hit in hits:
        sstart, send = hit["sstart"], hit["send"]
        qstart, qend = hit["qstart"], hit["qend"]
        slen = hit["slen"]

        start, end = min(sstart, send), max(sstart, send)
        strand = 1 if sstart <= send else -1

        # Expand query coords so subject spans 1..slen
        if start != 1:
            qstart -= (start - 1)
        if end != slen:
            qend += (slen - end)

        qstart0 = max(0, qstart - 1)
        contig_len = len(genome[hit["qaccver"]])
        qend_incl = min(qend, contig_len)

        if qend_incl <= qstart0:
            msg = (
                f"Invalid expanded coordinates for hit on contig {hit['qaccver']} "
                f"(qstart0={qstart0}, qend_incl={qend_incl}); skipping."
            )
            warnings.warn(msg, RuntimeWarning)
            print(f"[WARNING] {msg}", file=sys.stderr)
            continue

        alleles.append((hit["qaccver"], qstart0, qend_incl, hit["saccver"], strand))

    return alleles


def run_blast_locus(
    assembly_path: Path,
    locus_path: Path,
    genome: Fasta,
    assembly_name: str,
    min_cov_ratio: float,
) -> list:
    locus_name = locus_path.name.replace(".fasta", "")
    blast_output = os.path.join(LOCAL_CWD, f"blast_{assembly_name}_{locus_name}.txt")

    cmd = [
        "blastn",
        "-query", str(assembly_path),
        "-subject", str(locus_path),
        "-out", blast_output,
        "-outfmt", "6 qaccver saccver slen pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        "-max_target_seqs", "1000",
        "-num_threads", "1",
    ]

    start = datetime.now(timezone.utc)
    print(f"[{start.isoformat()}] Starting BLAST for locus={locus_name} assembly={assembly_name}")

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        end = datetime.now(timezone.utc)
        print(f"[{end.isoformat()}] BLAST FAILED for locus={locus_name} assembly={assembly_name} (elapsed: {end - start})")
        raise RuntimeError(
            f"BLAST failed for locus {locus_name} on {assembly_name}: {e}"
        ) from e
    else:
        end = datetime.now(timezone.utc)
        print(f"[{end.isoformat()}] Finished BLAST for locus={locus_name} assembly={assembly_name} (elapsed: {end - start})")

    alleles = []
    if os.path.exists(blast_output):
        alleles = parse_blast_output(blast_output, genome, min_cov_ratio)
        try:
            os.remove(blast_output)
        except OSError as e:
            msg = f"Could not remove temporary BLAST output {blast_output}: {e}"
            warnings.warn(msg, RuntimeWarning)
            print(f"[WARNING] {msg}", file=sys.stderr)

    return alleles


def process_single_assembly(
    assembly_path: Path,
    schema_dir: Path,
    output_file: Path,
    log: object,
    max_workers: int,
    # optional knobs (kept optional to avoid changing callers)
    max_stop_extend: int = 90,
    min_cov_ratio: float = 0.70,
) -> None:
    """
    New implementation that mirrors blast_locus.py behavior but keeps the original API.
    Writes a single FASTA to output_file.
    """
    if not assembly_path.exists():
        raise FileNotFoundError(f"Assembly file not found: {assembly_path}")

    try:
        genome = Fasta(str(assembly_path))
    except Exception as e:
        raise RuntimeError(f"Failed to open assembly FASTA with pyfaidx: {assembly_path}: {e}") from e

    loci = [f for f in os.listdir(schema_dir) if f.endswith(".fasta")]
    if not loci:
        raise ValueError(f"No .fasta locus files found in scheme directory: {schema_dir}")

    all_alleles = []
    assembly_name = assembly_path.stem

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                run_blast_locus,
                assembly_path,
                schema_dir / locus,
                genome,
                assembly_name,
                min_cov_ratio,
            ): locus for locus in loci
        }

        for future in as_completed(futures):
            locus = futures[future]
            try:
                all_alleles.extend(future.result())
            except Exception as e:
                # By default: fatal (matches your new script)
                raise RuntimeError(f"Failed processing locus {locus} for assembly {assembly_path.name}: {e}") from e

    extracted = extract_subsequences(genome, all_alleles, max_stop_extend)

    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w", encoding="utf-8") as out_f:
        for hdr, seq in extracted:
            out_f.write(f"{hdr}\n{seq}\n")

    print(f"Wrote {len(extracted)} alleles to {output_file}", file=sys.stderr)
