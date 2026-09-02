#!/usr/bin/env python
"""
Run per-locus BLAST searches against genome assemblies and extract valid CDSs.

Key behavior in this version:
- BLAST hits are grouped by physical genomic region.
- One representative allele hit is selected per qualifying region.
- Every distinct qualifying region is considered, including multiple/ambiguous
  copies for the same locus.
- Repeated matches from many scheme alleles to the same genomic coordinates
  produce only one extraction candidate.
- Candidates are compared with every allele in the locus in the same codon frame.
- A candidate that is a subsection of a longer known allele is never emitted.
  The script first tries to recover that complete allele from real assembly
  flanks; if exact recovery fails, the candidate is recorded in the stats TSV
  and excluded from the CDS FASTA.
- A longer candidate containing a shorter known allele is trimmed to the exact
  shorter allele sequence as it occurs in the assembly.
- A sequence differing from a known allele only at a valid start codon is
  normalized to the scheme allele for chewBBACA, with the observed and scheme
  start codons recorded in the header and stats TSV.
- No generic one- or two-base frame padding is performed.
- Stop-codon extension is disabled by default and is opt-in.
- Partial gapped hits are not used for speculative boundary extrapolation.
- Every candidate decision is written to a per-assembly TSV sidecar.
- Scheme FASTAs are parsed without shared .fai files.
"""

import csv
import logging
import os
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from shutil import copytree, ignore_patterns

from Bio import SeqIO
from Bio.Seq import Seq
from pyfaidx import Fasta


# Pick a base tmp dir that is *meant* for the current job/user
TMP_BASE = os.environ.get("TMPDIR", "/tmp")

# Create a unique, private temp dir for THIS run of the script
LOCAL_CWD = tempfile.mkdtemp(prefix="chewbbaca_runs_", dir=TMP_BASE)

# (optional) make sure permissions are sane even with umask weirdness
os.chmod(LOCAL_CWD, 0o700)

# Define valid bases and codon sets (uppercase)
VALID_BASES = {"A", "T", "C", "G"}
START_CODONS = {"ATG", "TTG", "GTG", "CTG", "ATA", "ATT"}
STOP_CODONS = {"TAA", "TAG", "TGA"}

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
    handlers=[logging.StreamHandler(sys.stderr)],
)



def is_valid_cds(seq: str) -> bool:
    """Return True when seq has a valid start, terminal stop, and codon length."""
    return (
        len(seq) >= 6
        and seq[:3] in START_CODONS
        and seq[-3:] in STOP_CODONS
        and len(seq) % 3 == 0
    )


def has_internal_stop(seq: str) -> bool:
    """Return True if seq contains an internal in-frame stop codon."""
    for i in range(3, len(seq) - 3, 3):
        if seq[i : i + 3] in STOP_CODONS:
            return True
    return False


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return str(Seq(seq).reverse_complement())


def only_start_codon_diff(query_seq: str, allele_seq: str) -> bool:
    if not query_seq or not allele_seq:
        return False

    query_seq = query_seq.upper()
    allele_seq = allele_seq.upper()

    return (
        len(query_seq) == len(allele_seq)
        and query_seq[3:] == allele_seq[3:]
        and query_seq[:3] in START_CODONS
        and allele_seq[:3] in START_CODONS
        and query_seq[:3] != allele_seq[:3]
    )


def allele_containment_case(query_seq: str, allele_seq: str) -> str:
    if not query_seq or not allele_seq:
        return ""

    query_seq = query_seq.upper()
    allele_seq = allele_seq.upper()

    if query_seq == allele_seq:
        return "exact"
    if query_seq in allele_seq:
        return "query_is_subsequence"
    if allele_seq in query_seq:
        return "allele_is_subsequence"
    return ""


def start_region_subsequence_case(query_seq: str, allele_seq: str) -> bool:
    """
    Return True only when query_seq is a suffix of allele_seq and exactly one
    valid start codon is missing from the query.
    """
    if not query_seq or not allele_seq:
        return False

    query_seq = query_seq.upper()
    allele_seq = allele_seq.upper()

    if query_seq == allele_seq or not allele_seq.endswith(query_seq):
        return False

    missing = len(allele_seq) - len(query_seq)
    return missing == 3 and allele_seq[:3] in START_CODONS


def allele_is_start_region_extension(query_seq: str, allele_seq: str) -> bool:
    if not query_seq or not allele_seq:
        return False

    query_seq = query_seq.upper()
    allele_seq = allele_seq.upper()

    if not query_seq.endswith(allele_seq):
        return False

    extra = len(query_seq) - len(allele_seq)
    return (
        extra == 3
        and query_seq[:3] in START_CODONS
        and allele_seq[:3] in START_CODONS
    )


def orient_sequence(
    seq: str,
    qstart: int,
    qend: int,
    strand: int,
):
    """Orient an assembly interval to the biological CDS direction."""
    stats = Counter()

    if strand == -1:
        seq = reverse_complement(seq)
        stats["reoriented"] += 1

    return seq.upper(), qstart, qend, stats


def validate_or_extend_novel_sequence(
    seq: str,
    qstart: int,
    qend: int,
    genome: Fasta,
    seq_id: str,
    header_raw: str,
    strand: int,
    max_stop_extend: int,
    alignment_gapped: bool,
):
    """
    Validate a sequence only after known-allele relationships were checked.

    No one- or two-base frame padding is performed. Optional stop extension is
    allowed only in whole-codon steps and only for an ungapped representative
    alignment.
    """
    stats = Counter()

    if len(seq) % 3 != 0:
        msg = (
            f"{header_raw}: extracted length {len(seq)} is not divisible by 3; "
            "skipping because frame padding is disabled."
        )
        print(f"[WARNING] {msg}", file=sys.stderr)
        logging.warning(msg)
        return (
            "", qstart, qend, stats,
            "excluded-out-of-frame",
            msg,
        )

    if seq[-3:] in STOP_CODONS:
        return seq, qstart, qend, stats, None, None

    if max_stop_extend == 0:
        msg = (
            f"{header_raw}: no terminal stop codon and stop extension is disabled; "
            "skipping."
        )
        print(f"[WARNING] {msg}", file=sys.stderr)
        logging.warning(msg)
        return (
            "", qstart, qend, stats,
            "excluded-no-terminal-stop",
            msg,
        )

    if alignment_gapped:
        msg = (
            f"{header_raw}: no terminal stop codon and the representative BLAST "
            "alignment is gapped; refusing speculative stop extension."
        )
        print(f"[WARNING] {msg}", file=sys.stderr)
        logging.warning(msg)
        return (
            "", qstart, qend, stats,
            "excluded-gapped-stop-extension",
            msg,
        )

    contig_len = len(genome[seq_id])

    for extension in range(3, max_stop_extend + 1, 3):
        if strand == 1:
            if qend + extension > contig_len:
                break
            cand = genome[seq_id][qstart : qend + extension].seq.upper()
        else:
            if qstart - extension < 0:
                break
            raw = genome[seq_id][qstart - extension : qend].seq.upper()
            cand = reverse_complement(raw)

        if len(cand) % 3 != 0:
            continue
        if cand[-3:] not in STOP_CODONS:
            continue
        if has_internal_stop(cand):
            continue

        if strand == 1:
            qend += extension
        else:
            qstart -= extension
        stats["stop_extended"] = extension
        return cand, qstart, qend, stats, None, None

    msg = (
        f"{header_raw}: no acceptable terminal stop within {max_stop_extend} bp "
        "of the biological 3' end; skipping."
    )
    print(f"[WARNING] {msg}", file=sys.stderr)
    logging.warning(msg)
    return (
        "", qstart, qend, stats,
        "excluded-no-acceptable-stop",
        msg,
    )

def hit_query_interval(rec: dict) -> tuple[int, int]:
    """Return the one-based inclusive genomic interval of a BLAST hit."""
    return min(rec["qstart"], rec["qend"]), max(rec["qstart"], rec["qend"])


def hit_strand(rec: dict) -> int:
    """
    Return the orientation of the allele relative to the forward query contig.

    A hit is forward when query and subject coordinates run in the same
    direction and reverse when they run in opposite directions.
    """
    query_direction = 1 if rec["qend"] >= rec["qstart"] else -1
    subject_direction = 1 if rec["send"] >= rec["sstart"] else -1
    return 1 if query_direction == subject_direction else -1


def subject_coverage(rec: dict) -> float:
    """Return the fraction of the subject allele covered by the alignment."""
    subject_span = abs(rec["send"] - rec["sstart"]) + 1
    return subject_span / rec["slen"] if rec["slen"] else 0.0


def is_exact_full_length(rec: dict) -> bool:
    """Return True for a clean, full-length, ungapped exact allele match."""
    query_span = abs(rec["qend"] - rec["qstart"]) + 1
    subject_span = abs(rec["send"] - rec["sstart"]) + 1
    return (
        rec["pident"] == 100.0
        and rec["mismatch"] == 0
        and rec["gapopen"] == 0
        and rec["length"] == rec["slen"]
        and query_span == rec["slen"]
        and subject_span == rec["slen"]
    )


def intervals_same_region(rec_a: dict, rec_b: dict, minimum_overlap: float) -> bool:
    """
    Return True when two hits most likely describe the same physical region.

    Hits must be on the same contig and strand. Their overlap must cover at
    least minimum_overlap of the shorter query interval. This groups shorter
    or longer allele variants at one coordinate while keeping separate copies
    on different contigs or non-overlapping positions distinct.
    """
    if rec_a["qaccver"] != rec_b["qaccver"]:
        return False
    if hit_strand(rec_a) != hit_strand(rec_b):
        return False

    a_start, a_end = hit_query_interval(rec_a)
    b_start, b_end = hit_query_interval(rec_b)

    overlap = max(0, min(a_end, b_end) - max(a_start, b_start) + 1)
    shorter = min(a_end - a_start + 1, b_end - b_start + 1)

    return shorter > 0 and overlap / shorter >= minimum_overlap


def cluster_blast_hits(records: list[dict], minimum_overlap: float) -> list[list[dict]]:
    """Group BLAST records into distinct physical genomic-region clusters."""
    clusters: list[list[dict]] = []

    sorted_records = sorted(
        records,
        key=lambda rec: (
            rec["qaccver"],
            hit_strand(rec),
            hit_query_interval(rec)[0],
            hit_query_interval(rec)[1],
        ),
    )

    for rec in sorted_records:
        matching_cluster = None

        for cluster in clusters:
            if any(
                intervals_same_region(rec, member, minimum_overlap)
                for member in cluster
            ):
                matching_cluster = cluster
                break

        if matching_cluster is None:
            clusters.append([rec])
        else:
            matching_cluster.append(rec)

    return clusters


def representative_hit_key(rec: dict) -> tuple:
    """
    Rank allele hits only within one genomic region.

    A full-subject or ungapped hit is preferred when boundary inference may be
    required. This is not a blanket rule that every ungapped hit beats every
    gapped hit: full-length gapped hits can still rank by identity.
    """
    normalized_bitscore = rec["bitscore"] / max(rec["length"], 1)
    coverage = subject_coverage(rec)
    boundary_reliable = coverage == 1.0 or rec["gapopen"] == 0

    return (
        is_exact_full_length(rec),
        boundary_reliable,
        coverage,
        rec["pident"],
        rec["gapopen"] == 0,
        -rec.get("gaps", 0),
        -rec["gapopen"],
        -rec["mismatch"],
        normalized_bitscore,
        rec["bitscore"],
    )


def select_candidate_regions(
    records: list[dict],
    min_cov_ratio: float,
    min_identity: float,
    region_overlap: float,
    locus_name: str,
) -> list[dict]:
    """
    Select one representative allele hit for every qualifying genomic region.

    No global single-best fallback is used. Therefore, multiple physical copies
    are retained, but multiple allele references matching one physical copy are
    reduced to one representative.
    """
    clusters = cluster_blast_hits(records, region_overlap)
    candidates = []

    for cluster in clusters:
        qualifying_hits = [
            rec
            for rec in cluster
            if subject_coverage(rec) >= min_cov_ratio
            and rec["pident"] >= min_identity
        ]

        if not qualifying_hits:
            continue

        representative = max(qualifying_hits, key=representative_hit_key).copy()
        representative["cluster_hit_count"] = len(cluster)
        representative["qualifying_hit_count"] = len(qualifying_hits)
        candidates.append(representative)

    candidates.sort(
        key=lambda rec: (
            rec["qaccver"],
            hit_query_interval(rec)[0],
            hit_query_interval(rec)[1],
            rec["saccver"],
        )
    )

    candidate_count = len(candidates)
    for candidate_number, rec in enumerate(candidates, start=1):
        rec["candidate_number"] = candidate_number
        rec["candidate_count"] = candidate_count

    if candidate_count == 0:
        logging.info(
            "%s: no genomic region passed min_cov=%.3f and min_identity=%.3f.",
            locus_name,
            min_cov_ratio,
            min_identity,
        )
    elif candidate_count == 1:
        rec = candidates[0]
        logging.info(
            "%s: retained one genomic candidate at %s:%d-%d; nearest allele=%s, "
            "identity=%.3f, subject_coverage=%.3f, exact=%s.",
            locus_name,
            rec["qaccver"],
            hit_query_interval(rec)[0],
            hit_query_interval(rec)[1],
            rec["saccver"],
            rec["pident"],
            subject_coverage(rec),
            is_exact_full_length(rec),
        )
    else:
        summary = "; ".join(
            (
                f"candidate-{rec['candidate_number']}="
                f"{rec['qaccver']}:{hit_query_interval(rec)[0]}-"
                f"{hit_query_interval(rec)[1]},nearest={rec['saccver']},"
                f"identity={rec['pident']:.3f},coverage={subject_coverage(rec):.3f},"
                f"exact={is_exact_full_length(rec)}"
            )
            for rec in candidates
        )
        logging.warning(
            "%s: retained %d distinct genomic candidates for multi-copy/ambiguous "
            "handling: %s",
            locus_name,
            candidate_count,
            summary,
        )

    return candidates


def expand_hit_to_subject_boundaries(
    hit: dict,
    contig_len: int,
):
    """
    Expand query coordinates to the expected subject-allele boundaries.

    Returns:
      qstart0, qend_exclusive, left_truncated, right_truncated,
      used_boundary_extrapolation, unsafe_gapped_extrapolation

    Expected coordinates outside the contig are reported rather than silently
    clamped. Boundary extrapolation from a partial gapped alignment is marked
    unsafe because subject and query offsets are not reliably one-to-one.
    """
    q_low, q_high = hit_query_interval(hit)
    s_low = min(hit["sstart"], hit["send"])
    s_high = max(hit["sstart"], hit["send"])

    missing_subject_start = s_low - 1
    missing_subject_end = hit["slen"] - s_high
    used_boundary_extrapolation = bool(
        missing_subject_start or missing_subject_end
    )

    unsafe_gapped_extrapolation = (
        used_boundary_extrapolation and hit["gapopen"] > 0
    )

    if hit_strand(hit) == 1:
        requested_low = q_low - missing_subject_start
        requested_high = q_high + missing_subject_end
    else:
        requested_low = q_low - missing_subject_end
        requested_high = q_high + missing_subject_start

    left_truncated = requested_low < 1
    right_truncated = requested_high > contig_len

    qstart0 = max(1, requested_low) - 1
    qend_exclusive = min(contig_len, requested_high)

    return (
        qstart0,
        qend_exclusive,
        left_truncated,
        right_truncated,
        used_boundary_extrapolation,
        unsafe_gapped_extrapolation,
    )

def parse_blast_output(
    blast_output: str,
    genome: Fasta,
    min_cov_ratio: float,
    min_identity: float,
    region_overlap: float,
    locus_name: str,
    locus_path: str,
) -> tuple[list, list[dict]]:
    """
    Parse BLAST output and return extraction records plus early decisions.

    Early decisions cover loci with no qualifying region and candidates that
    cannot safely proceed to sequence extraction because expected boundaries
    cross a contig edge or require gapped coordinate extrapolation.
    """
    records = []
    decisions = []

    with open(blast_output) as fh:
        for line_number, line in enumerate(fh, start=1):
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 14:
                logging.warning(
                    "%s line %d: expected 14 BLAST columns, found %d; skipping.",
                    blast_output,
                    line_number,
                    len(cols),
                )
                continue

            try:
                rec = {
                    "qaccver": cols[0],
                    "saccver": cols[1],
                    "slen": int(cols[2]),
                    "pident": float(cols[3]),
                    "length": int(cols[4]),
                    "mismatch": int(cols[5]),
                    "gapopen": int(cols[6]),
                    "gaps": int(cols[7]),
                    "qstart": int(cols[8]),
                    "qend": int(cols[9]),
                    "sstart": int(cols[10]),
                    "send": int(cols[11]),
                    "evalue": float(cols[12]),
                    "bitscore": float(cols[13]),
                }
            except ValueError as exc:
                logging.warning(
                    "%s line %d: invalid numeric BLAST value (%s); skipping.",
                    blast_output,
                    line_number,
                    exc,
                )
                continue

            if rec["qaccver"] not in genome:
                logging.warning(
                    "%s: query contig %s is absent from the assembly FASTA; skipping hit.",
                    locus_name,
                    rec["qaccver"],
                )
                continue

            records.append(rec)

    hits = select_candidate_regions(
        records=records,
        min_cov_ratio=min_cov_ratio,
        min_identity=min_identity,
        region_overlap=region_overlap,
        locus_name=locus_name,
    )

    if not hits:
        decisions.append(
            make_decision_row(
                locus_name=locus_name,
                status="not-found-no-qualifying-region",
                reason=(
                    f"No BLAST region passed min_cov={min_cov_ratio:.3f} "
                    f"and min_identity={min_identity:.3f}."
                ),
            )
        )
        return [], decisions

    alleles = []

    for hit in hits:
        contig_len = len(genome[hit["qaccver"]])
        (
            qstart0,
            qend_exclusive,
            left_truncated,
            right_truncated,
            used_boundary_extrapolation,
            unsafe_gapped_extrapolation,
        ) = expand_hit_to_subject_boundaries(hit, contig_len)
        strand = hit_strand(hit)
        common = dict(
            locus_name=locus_name,
            candidate_number=hit["candidate_number"],
            candidate_count=hit["candidate_count"],
            seq_id=hit["qaccver"],
            start=max(1, qstart0 + 1),
            end=min(contig_len, qend_exclusive),
            strand=strand,
            nearest_allele=hit["saccver"],
            pident=hit["pident"],
            coverage=subject_coverage(hit),
            blast_exact=is_exact_full_length(hit),
            gapopen=hit["gapopen"],
            gaps=hit.get("gaps", 0),
            boundary_extrapolated=used_boundary_extrapolation,
        )

        if left_truncated or right_truncated:
            reason = (
                "Expected allele boundaries extend outside the contig; "
                f"left_truncated={left_truncated}, "
                f"right_truncated={right_truncated}. The remainder is unknown."
            )
            logging.warning(
                "%s candidate %d: %s",
                locus_name,
                hit["candidate_number"],
                reason,
            )
            decisions.append(
                make_decision_row(
                    **common,
                    status="not-found-partial-contig-edge",
                    reason=reason,
                )
            )
            continue

        if unsafe_gapped_extrapolation:
            reason = (
                "Partial representative alignment is gapped, so missing allele "
                "boundaries cannot be inferred reliably."
            )
            logging.warning(
                "%s candidate %d: %s",
                locus_name,
                hit["candidate_number"],
                reason,
            )
            decisions.append(
                make_decision_row(
                    **common,
                    status="excluded-gapped-boundary-extrapolation",
                    reason=reason,
                )
            )
            continue

        if qend_exclusive <= qstart0:
            reason = "Expanded genomic interval is empty or reversed."
            logging.warning(
                "%s candidate %d: %s",
                locus_name,
                hit["candidate_number"],
                reason,
            )
            decisions.append(
                make_decision_row(
                    **common,
                    status="excluded-invalid-interval",
                    reason=reason,
                )
            )
            continue

        alleles.append(
            (
                hit["qaccver"],
                qstart0,
                qend_exclusive,
                hit["saccver"],
                strand,
                locus_name,
                locus_path,
                hit["candidate_number"],
                hit["candidate_count"],
                hit["pident"],
                subject_coverage(hit),
                is_exact_full_length(hit),
                hit["gapopen"],
                hit.get("gaps", 0),
                used_boundary_extrapolation,
            )
        )

    return alleles, decisions


def build_locus_allele_catalog(locus_path: str) -> tuple:
    """
    Load all alleles for one locus without creating or reading a shared .fai.

    Scheme locus files are small enough to parse sequentially with Biopython.
    Avoiding pyfaidx here prevents concurrent samples from racing while they
    create the same locus .fai file, which can leave a truncated index and
    abort another sample.

    Returns:
      allele_entries: list of (allele_id, sequence)
      allele_lookup: allele_id -> sequence
      exact_index: complete sequence -> sorted allele IDs
      coding_tail_index: (length, sequence after first codon) -> allele IDs
    """
    allele_entries = []
    allele_lookup = {}
    exact_index = defaultdict(list)
    coding_tail_index = defaultdict(list)

    for record in SeqIO.parse(locus_path, "fasta"):
        allele_id = record.id
        sequence = str(record.seq).upper()

        if allele_id in allele_lookup:
            raise ValueError(
                f"Duplicate allele ID {allele_id!r} in locus FASTA {locus_path}"
            )
        if not sequence:
            raise ValueError(
                f"Empty allele sequence for {allele_id!r} in {locus_path}"
            )

        allele_entries.append((allele_id, sequence))
        allele_lookup[allele_id] = sequence
        exact_index[sequence].append(allele_id)

        if len(sequence) >= 3:
            coding_tail_index[(len(sequence), sequence[3:])].append(allele_id)

    if not allele_entries:
        raise ValueError(f"No FASTA records found in locus file: {locus_path}")

    allele_entries.sort(key=lambda item: item[0])
    for allele_ids in exact_index.values():
        allele_ids.sort()
    for allele_ids in coding_tail_index.values():
        allele_ids.sort()

    return allele_entries, allele_lookup, exact_index, coding_tail_index


def prefer_allele_id(allele_ids: list[str], nearest_saccver: str) -> str:
    """Choose the nearest BLAST allele when possible, otherwise deterministically."""
    if nearest_saccver in allele_ids:
        return nearest_saccver
    return sorted(allele_ids)[0]


def find_start_codon_difference(
    query_seq: str,
    allele_lookup: dict,
    coding_tail_index: dict,
    nearest_saccver: str,
):
    """Return an allele differing only in its valid start codon, if present."""
    if len(query_seq) < 3:
        return None

    allele_ids = coding_tail_index.get((len(query_seq), query_seq[3:]), [])
    allele_ids = [
        allele_id
        for allele_id in allele_ids
        if query_seq[:3] in START_CODONS
        and allele_lookup[allele_id][:3] in START_CODONS
        and query_seq[:3] != allele_lookup[allele_id][:3]
    ]

    if not allele_ids:
        return None

    allele_id = prefer_allele_id(allele_ids, nearest_saccver)
    return allele_id, allele_lookup[allele_id]


def find_in_frame_occurrences(shorter: str, longer: str) -> list[int]:
    """Return all zero-based occurrences of shorter in longer in codon frame."""
    positions = []
    position = longer.find(shorter)

    while position != -1:
        if position % 3 == 0:
            positions.append(position)
        position = longer.find(shorter, position + 1)

    return positions


def find_containing_alleles(
    query_seq: str,
    allele_entries: list[tuple[str, str]],
    nearest_saccver: str,
) -> list[tuple[str, str, int, int]]:
    """
    Return longer known alleles containing query_seq in the same codon frame.

    Results are ordered by the smallest amount of missing sequence, then by the
    nearest BLAST allele and allele ID. Each tuple contains:
      allele_id, allele_sequence, missing_bp, query_offset_in_allele
    """
    matches = []

    for allele_id, allele_seq in allele_entries:
        if len(allele_seq) <= len(query_seq):
            continue

        positions = find_in_frame_occurrences(query_seq, allele_seq)
        for position in positions:
            missing = len(allele_seq) - len(query_seq)
            nearest_penalty = 0 if allele_id == nearest_saccver else 1
            matches.append(
                (
                    missing,
                    nearest_penalty,
                    allele_id,
                    position,
                    allele_seq,
                )
            )

    matches.sort()
    return [
        (allele_id, allele_seq, missing, position)
        for missing, _, allele_id, position, allele_seq in matches
    ]


def recover_known_allele_from_assembly(
    genome: Fasta,
    seq_id: str,
    qstart: int,
    qend: int,
    strand: int,
    query_seq: str,
    allele_seq: str,
    query_offset_in_allele: int,
):
    """
    Try to recover a complete known allele using only real assembly flanks.

    qstart/qend describe query_seq in zero-based, end-exclusive genomic
    coordinates. query_offset_in_allele is the in-frame position of query_seq
    inside allele_seq.

    Returns a dictionary with:
      recovered, sequence, start, end, touches_contig_edge, reason
    """
    missing_before = query_offset_in_allele
    missing_after = len(allele_seq) - (
        query_offset_in_allele + len(query_seq)
    )

    if missing_before < 0 or missing_after < 0:
        raise ValueError("Invalid containment offsets during allele recovery")

    if strand == 1:
        requested_start = qstart - missing_before
        requested_end = qend + missing_after
    else:
        requested_start = qstart - missing_after
        requested_end = qend + missing_before

    contig_len = len(genome[seq_id])
    touches_contig_edge = (
        requested_start < 0 or requested_end > contig_len
    )

    if touches_contig_edge:
        return {
            "recovered": False,
            "sequence": None,
            "start": requested_start,
            "end": requested_end,
            "touches_contig_edge": True,
            "reason": (
                "The complete known allele would extend outside the assembly "
                "contig."
            ),
        }

    raw = genome[seq_id][requested_start:requested_end].seq.upper()
    invalid = [base for base in raw if base not in VALID_BASES]
    if invalid:
        return {
            "recovered": False,
            "sequence": None,
            "start": requested_start,
            "end": requested_end,
            "touches_contig_edge": False,
            "reason": "The required assembly flanks contain ambiguous bases.",
        }

    oriented = raw if strand == 1 else reverse_complement(raw)
    if oriented != allele_seq:
        return {
            "recovered": False,
            "sequence": oriented,
            "start": requested_start,
            "end": requested_end,
            "touches_contig_edge": False,
            "reason": (
                "Assembly flanks are present but do not exactly recreate the "
                "known longer allele."
            ),
        }

    return {
        "recovered": True,
        "sequence": oriented,
        "start": requested_start,
        "end": requested_end,
        "touches_contig_edge": False,
        "reason": "Recovered exact known allele from assembly flanks.",
    }

def find_contained_allele(
    query_seq: str,
    allele_entries: list[tuple[str, str]],
    nearest_saccver: str,
):
    """
    Find the longest shorter known allele contained in query in the same frame.

    The longest match is preferred to minimize trimming. The nearest BLAST
    allele and then allele ID are deterministic tie-breakers.
    """
    matches = []

    for allele_id, allele_seq in allele_entries:
        if len(allele_seq) >= len(query_seq):
            continue

        positions = find_in_frame_occurrences(allele_seq, query_seq)
        for position in positions:
            nearest_penalty = 0 if allele_id == nearest_saccver else 1
            matches.append(
                (
                    -len(allele_seq),
                    nearest_penalty,
                    allele_id,
                    position,
                    allele_seq,
                )
            )

    if not matches:
        return None

    _, _, allele_id, position, allele_seq = min(matches)
    return allele_id, allele_seq, position


def trim_oriented_interval_to_subsequence(
    qstart: int,
    qend: int,
    strand: int,
    offset: int,
    subseq_length: int,
) -> tuple[int, int]:
    """Map an oriented-sequence substring back to genomic slice coordinates."""
    if strand == 1:
        new_start = qstart + offset
        new_end = new_start + subseq_length
    else:
        new_start = qend - (offset + subseq_length)
        new_end = qend - offset

    return new_start, new_end


def fasta_token(value: str) -> str:
    """Make a value safe for the first whitespace-delimited FASTA token."""
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_") or "unnamed"

STATS_FIELDS = [
    "assembly",
    "locus",
    "candidate",
    "candidate_number",
    "candidate_count",
    "included_in_cds_fasta",
    "output_id",
    "status",
    "reason",
    "contig",
    "start",
    "end",
    "strand",
    "nearest_allele",
    "matched_allele",
    "pident",
    "coverage",
    "blast_exact",
    "gapopen",
    "gaps",
    "boundary_extrapolated",
    "observed_start",
    "scheme_start",
    "output_length",
]


def make_decision_row(
    *,
    locus_name: str,
    candidate_number=None,
    candidate_count=None,
    included=False,
    output_id="",
    status: str,
    reason: str,
    seq_id="",
    start=None,
    end=None,
    strand=None,
    nearest_allele="",
    matched_allele="",
    pident=None,
    coverage=None,
    blast_exact=None,
    gapopen=None,
    gaps=None,
    boundary_extrapolated=None,
    observed_start="",
    scheme_start="",
    output_length=None,
) -> dict:
    """Create one normalized candidate-decision record for the stats TSV."""
    candidate = (
        f"candidate-{candidate_number}-of-{candidate_count}"
        if candidate_number is not None and candidate_count is not None
        else ""
    )

    return {
        "locus": locus_name,
        "candidate": candidate,
        "candidate_number": candidate_number,
        "candidate_count": candidate_count,
        "included_in_cds_fasta": "yes" if included else "no",
        "output_id": output_id,
        "status": status,
        "reason": reason,
        "contig": seq_id,
        "start": start,
        "end": end,
        "strand": (
            "+" if strand == 1 else "-" if strand == -1 else ""
        ),
        "nearest_allele": nearest_allele,
        "matched_allele": matched_allele,
        "pident": pident,
        "coverage": coverage,
        "blast_exact": blast_exact,
        "gapopen": gapopen,
        "gaps": gaps,
        "boundary_extrapolated": boundary_extrapolated,
        "observed_start": observed_start,
        "scheme_start": scheme_start,
        "output_length": output_length,
    }


def write_decision_stats(
    stats_path: str,
    assembly_name: str,
    decisions: list[dict],
) -> None:
    """Write all locus/candidate decisions for one assembly to a TSV file."""
    rows = sorted(
        decisions,
        key=lambda row: (
            row.get("locus", ""),
            row.get("candidate_number")
            if row.get("candidate_number") is not None
            else 0,
            row.get("status", ""),
        ),
    )

    with open(stats_path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=STATS_FIELDS,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()

        for row in rows:
            complete = {field: "" for field in STATS_FIELDS}
            complete.update(row)
            complete["assembly"] = assembly_name

            for key, value in list(complete.items()):
                if value is None:
                    complete[key] = ""

            writer.writerow(complete)


def extract_subsequences(
    genome: Fasta,
    alleles: list,
    max_stop_extend: int,
    assembly_name: str,
) -> tuple[list, list[dict]]:
    """
    Validate and extract every distinct locus/genomic-region candidate.

    Policy:
    - exact known allele: emit the assembly sequence;
    - valid start-codon-only difference: emit the scheme allele sequence for
      normalization, while recording observed and scheme start codons;
    - candidate contained in a longer known allele: try exact recovery from
      real assembly flanks; if recovery fails, never emit the subsection;
    - candidate containing a shorter known allele: trim to the longest exact
      in-frame known allele as it occurs in the assembly;
    - otherwise validate and emit the assembly sequence as novel.
    """
    extracted = []
    decisions = []
    seen_regions = set()
    alleles_by_locus = defaultdict(list)

    for allele in alleles:
        locus_path = allele[6]
        alleles_by_locus[locus_path].append(allele)

    assembly_token = fasta_token(assembly_name)

    for locus_path in sorted(alleles_by_locus):
        try:
            (
                allele_entries,
                allele_lookup,
                exact_index,
                coding_tail_index,
            ) = build_locus_allele_catalog(locus_path)

            locus_candidates = sorted(
                alleles_by_locus[locus_path],
                key=lambda item: (
                    item[5],
                    item[0],
                    item[1],
                    item[2],
                    item[7],
                ),
            )

            for allele in locus_candidates:
                (
                    seq_id,
                    qstart,
                    qend,
                    saccver,
                    strand,
                    locus_name,
                    _locus_path,
                    candidate_number,
                    candidate_count,
                    pident,
                    coverage,
                    blast_exact,
                    gapopen,
                    gaps,
                    used_boundary_extrapolation,
                ) = allele

                candidate_label = (
                    f"candidate-{candidate_number}-of-{candidate_count}"
                )
                # Keep the nearest/representative allele in the first
                # whitespace-delimited FASTA token. candidate_number preserves
                # uniqueness when multiple genomic copies have the same nearest
                # allele.
                record_id = (
                    f"{assembly_token}__{fasta_token(saccver)}__cand"
                    f"{candidate_number}"
                )
                header_raw = f">{record_id}"

                def add_decision(
                    status: str,
                    reason: str,
                    *,
                    included: bool = False,
                    output_id: str = "",
                    matched_allele: str = "",
                    start: int | None = None,
                    end: int | None = None,
                    observed_start: str = "",
                    scheme_start: str = "",
                    output_length: int | None = None,
                ) -> None:
                    decisions.append(
                        make_decision_row(
                            locus_name=locus_name,
                            candidate_number=candidate_number,
                            candidate_count=candidate_count,
                            included=included,
                            output_id=output_id,
                            status=status,
                            reason=reason,
                            seq_id=seq_id,
                            start=(qstart + 1 if start is None else start),
                            end=(qend if end is None else end),
                            strand=strand,
                            nearest_allele=saccver,
                            matched_allele=matched_allele,
                            pident=pident,
                            coverage=coverage,
                            blast_exact=blast_exact,
                            gapopen=gapopen,
                            gaps=gaps,
                            boundary_extrapolated=used_boundary_extrapolation,
                            observed_start=observed_start,
                            scheme_start=scheme_start,
                            output_length=output_length,
                        )
                    )

                subseq = genome[seq_id][qstart:qend].seq

                if len(subseq) != qend - qstart:
                    reason = (
                        f"Slice length {len(subseq)} does not equal expected "
                        f"length {qend - qstart}."
                    )
                    logging.error("%s: %s", header_raw, reason)
                    add_decision("excluded-slice-length-error", reason)
                    continue

                invalid_positions = [
                    (i, base)
                    for i, base in enumerate(subseq)
                    if base.upper() not in VALID_BASES
                ]
                if invalid_positions:
                    reason = (
                        f"Assembly interval contains {len(invalid_positions)} "
                        "ambiguous or invalid base(s)."
                    )
                    logging.warning("%s: %s", header_raw, reason)
                    add_decision("excluded-ambiguous-bases", reason)
                    continue

                (
                    fixed,
                    new_start,
                    new_end,
                    adjustment_stats,
                ) = orient_sequence(
                    subseq,
                    qstart,
                    qend,
                    strand,
                )

                output_seq = None
                output_status = None
                output_reason = None
                matched_saccver = None
                observed_start = ""
                scheme_start = ""
                sequence_source = "assembly"

                # 1) Exact sequence match against any known allele.
                exact_matches = exact_index.get(fixed, [])
                if exact_matches:
                    matched_saccver = prefer_allele_id(exact_matches, saccver)
                    output_seq = fixed
                    output_status = "exact"
                    output_reason = "Assembly sequence exactly matches a scheme allele."

                # 2) Normalize a valid start-codon-only difference to the scheme
                #    allele so chewBBACA does not create an artificial new allele.
                if output_seq is None:
                    start_match = find_start_codon_difference(
                        fixed,
                        allele_lookup,
                        coding_tail_index,
                        saccver,
                    )
                    if start_match is not None:
                        matched_saccver, matched_seq = start_match
                        observed_start = fixed[:3]
                        scheme_start = matched_seq[:3]
                        output_seq = matched_seq
                        output_status = "normalized-start-codon"
                        output_reason = (
                            f"Assembly differs from {matched_saccver} only at "
                            f"the valid start codon ({observed_start} -> "
                            f"{scheme_start}); emitted scheme allele sequence."
                        )
                        sequence_source = "scheme-start-normalization"
                        logging.info("%s: %s", header_raw, output_reason)

                # 3) A subsection of a longer known allele is never emitted.
                #    First try to recover an exact complete known allele using
                #    actual assembly flanks. No scheme bases are inserted.
                if output_seq is None:
                    containing_matches = find_containing_alleles(
                        fixed,
                        allele_entries,
                        saccver,
                    )
                    if containing_matches:
                        recovered_match = None
                        recovery_attempts = []

                        for (
                            allele_id,
                            allele_seq,
                            missing,
                            position,
                        ) in containing_matches:
                            recovery = recover_known_allele_from_assembly(
                                genome=genome,
                                seq_id=seq_id,
                                qstart=new_start,
                                qend=new_end,
                                strand=strand,
                                query_seq=fixed,
                                allele_seq=allele_seq,
                                query_offset_in_allele=position,
                            )
                            recovery_attempts.append(
                                (allele_id, missing, position, recovery)
                            )
                            if recovery["recovered"]:
                                recovered_match = (
                                    allele_id,
                                    allele_seq,
                                    missing,
                                    position,
                                    recovery,
                                )
                                break

                        if recovered_match is not None:
                            (
                                matched_saccver,
                                _matched_seq,
                                missing,
                                position,
                                recovery,
                            ) = recovered_match
                            output_seq = recovery["sequence"]
                            new_start = recovery["start"]
                            new_end = recovery["end"]
                            output_status = "recovered-known-allele"
                            output_reason = (
                                f"Candidate was an in-frame subsection of "
                                f"{matched_saccver}; recovered the complete exact "
                                f"allele from {missing} bp of real assembly flanks."
                            )
                            sequence_source = "assembly-flank-recovery"
                            logging.info("%s: %s", header_raw, output_reason)
                        else:
                            all_edge = all(
                                attempt[3]["touches_contig_edge"]
                                for attempt in recovery_attempts
                            )
                            closest_id, missing, position, closest = recovery_attempts[0]
                            location = (
                                "prefix"
                                if position == 0
                                else "suffix"
                                if position + len(fixed)
                                == len(allele_lookup[closest_id])
                                else "internal subsection"
                            )

                            if all_edge:
                                status = "not-found-partial-contig-edge"
                                reason = (
                                    f"Candidate is an in-frame {location} of "
                                    f"known allele {closest_id}, missing {missing} "
                                    "bp, but the complete allele would cross a "
                                    "contig edge. The unknown remainder is not "
                                    "added to the CDS FASTA."
                                )
                            else:
                                status = "excluded-unresolved-known-subsection"
                                reason = (
                                    f"Candidate is an in-frame {location} of "
                                    f"known allele {closest_id}, missing {missing} "
                                    "bp. Assembly flanks did not exactly recover "
                                    "a complete known allele, so the subsection is "
                                    "not emitted as a novel CDS."
                                )

                            logging.info("%s: %s", header_raw, reason)
                            add_decision(
                                status,
                                reason,
                                matched_allele=closest_id,
                                start=new_start + 1,
                                end=new_end,
                            )
                            continue

                # 4) A longer candidate containing a shorter exact scheme allele
                #    is trimmed to the longest known allele substring found in
                #    the assembly. This avoids calling uncertain extra sequence as
                #    a distinct allele.
                if output_seq is None:
                    contained_match = find_contained_allele(
                        fixed,
                        allele_entries,
                        saccver,
                    )
                    if contained_match is not None:
                        matched_saccver, matched_seq, offset = contained_match
                        trimmed_seq = fixed[offset : offset + len(matched_seq)]

                        if trimmed_seq != matched_seq:
                            raise AssertionError(
                                "Contained-allele lookup returned inconsistent "
                                "sequence coordinates."
                            )

                        new_start, new_end = trim_oriented_interval_to_subsequence(
                            new_start,
                            new_end,
                            strand,
                            offset,
                            len(trimmed_seq),
                        )
                        output_seq = trimmed_seq
                        output_status = "trimmed-to-known"
                        output_reason = (
                            f"Longer candidate contains exact in-frame allele "
                            f"{matched_saccver}; emitted that assembly substring "
                            "instead of creating a boundary-derived allele."
                        )
                        sequence_source = "assembly-trimmed"
                        logging.info("%s: %s", header_raw, output_reason)

                # 5) Validate an otherwise novel assembly-derived CDS only
                #    after exact, start-normalization, and containment checks.
                if output_seq is None:
                    (
                        fixed,
                        new_start,
                        new_end,
                        novel_adjustments,
                        failure_status,
                        failure_reason,
                    ) = validate_or_extend_novel_sequence(
                        fixed,
                        new_start,
                        new_end,
                        genome,
                        seq_id,
                        header_raw,
                        strand,
                        max_stop_extend,
                        alignment_gapped=gapopen > 0,
                    )
                    adjustment_stats.update(novel_adjustments)

                    if not fixed:
                        add_decision(
                            failure_status or "excluded-sequence-validation",
                            failure_reason or "Sequence validation failed.",
                            start=new_start + 1,
                            end=new_end,
                        )
                        continue

                    if has_internal_stop(fixed):
                        reason = (
                            "Assembly sequence contains an internal in-frame stop "
                            "codon."
                        )
                        logging.warning("%s: %s", header_raw, reason)
                        add_decision(
                            "excluded-internal-stop",
                            reason,
                            start=new_start + 1,
                            end=new_end,
                        )
                        continue

                    if not is_valid_cds(fixed):
                        reason = (
                            "Assembly sequence lacks an accepted start/stop or "
                            "has a non-codon length."
                        )
                        logging.warning("%s: %s", header_raw, reason)
                        add_decision(
                            "excluded-invalid-cds",
                            reason,
                            start=new_start + 1,
                            end=new_end,
                        )
                        continue

                    output_seq = fixed
                    output_status = "novel"
                    output_reason = (
                        "Valid assembly-derived CDS with no exact or containment "
                        "relationship to a known allele."
                    )
                    sequence_source = "assembly"

                if has_internal_stop(output_seq) or not is_valid_cds(output_seq):
                    reason = (
                        "Final selected sequence failed CDS validation after "
                        "normalization or containment handling."
                    )
                    logging.warning("%s: %s", header_raw, reason)
                    add_decision(
                        "excluded-final-cds-validation",
                        reason,
                        matched_allele=matched_saccver or "",
                        start=new_start + 1,
                        end=new_end,
                    )
                    continue

                region_key = (locus_name, seq_id, new_start, new_end, strand)
                if region_key in seen_regions:
                    reason = (
                        "Final genomic region duplicates an already emitted "
                        "candidate for this locus."
                    )
                    logging.warning("%s: %s", header_raw, reason)
                    add_decision(
                        "excluded-duplicate-region",
                        reason,
                        matched_allele=matched_saccver or "",
                        start=new_start + 1,
                        end=new_end,
                    )
                    continue
                seen_regions.add(region_key)

                description_parts = [
                    f"locus={locus_name}",
                    f"blast_candidate={candidate_number}/{candidate_count}",
                    f"status={output_status}",
                    f"source={sequence_source}",
                    f"contig={seq_id}",
                    f"coords={new_start + 1}-{new_end}",
                    f"strand={'+' if strand == 1 else '-'}",
                    f"nearest={saccver}",
                    f"matched={matched_saccver or 'none'}",
                    f"pident={pident:.3f}",
                    f"coverage={coverage:.3f}",
                    f"gapopen={gapopen}",
                    f"gaps={gaps}",
                    f"boundary_extrapolated={used_boundary_extrapolation}",
                ]
                if observed_start:
                    description_parts.append(f"observed_start={observed_start}")
                if scheme_start:
                    description_parts.append(f"scheme_start={scheme_start}")

                final_hdr = f">{record_id} " + " ".join(description_parts)
                extracted.append((final_hdr, output_seq))

                add_decision(
                    output_status,
                    output_reason or "Included in CDS FASTA.",
                    included=True,
                    output_id=record_id,
                    matched_allele=matched_saccver or "",
                    start=new_start + 1,
                    end=new_end,
                    observed_start=observed_start,
                    scheme_start=scheme_start,
                    output_length=len(output_seq),
                )

                logging.info(
                    "%s: extracted status=%s; matched_allele=%s, BLAST "
                    "identity=%.3f, subject_coverage=%.3f, blast_exact=%s, "
                    "final_length=%d, adjustments=%s.",
                    final_hdr,
                    output_status,
                    matched_saccver or "none",
                    pident,
                    coverage,
                    blast_exact,
                    len(output_seq),
                    dict(adjustment_stats),
                )
        except Exception:
            logging.exception(
                "Failed while reading or evaluating scheme locus %s.",
                locus_path,
            )
            raise

    return extracted, decisions


def run_blast_locus(
    assembly_path: str,
    locus_path: str,
    genome: Fasta,
    assembly_name: str,
    min_cov_ratio: float,
    min_identity: float,
    region_overlap: float,
) -> tuple[list, list[dict]]:
    """Run BLAST for one locus and return candidates plus early decisions."""
    locus_name = os.path.basename(locus_path).replace(".fasta", "")
    blast_output = os.path.join(
        LOCAL_CWD, f"blast_{assembly_name}_{locus_name}.txt"
    )

    cmd = [
        "blastn",
        "-query",
        assembly_path,
        "-subject",
        locus_path,
        "-out",
        blast_output,
        "-outfmt",
        (
            "6 qaccver saccver slen pident length mismatch gapopen gaps "
            "qstart qend sstart send evalue bitscore"
        ),
        "-max_target_seqs",
        "1000",
        "-num_threads",
        "1",
    ]

    blast_started_at = datetime.now(timezone.utc)
    blast_started_perf = time.perf_counter()
    logging.info(
        "[%s] Starting BLAST for locus=%s assembly=%s",
        blast_started_at.isoformat(),
        locus_name,
        assembly_name,
    )

    alleles = []
    decisions = []

    try:
        subprocess.run(cmd, check=True)

        # Parsing is intentionally included in the per-locus elapsed time.
        # This makes the timing representative of the complete BLAST-locus
        # operation used by the pipeline, not only the external blastn process.
        if os.path.exists(blast_output):
            alleles, decisions = parse_blast_output(
                blast_output=blast_output,
                genome=genome,
                min_cov_ratio=min_cov_ratio,
                min_identity=min_identity,
                region_overlap=region_overlap,
                locus_name=locus_name,
                locus_path=locus_path,
            )
        else:
            decisions.append(
                make_decision_row(
                    locus_name=locus_name,
                    status="not-found-no-blast-output",
                    reason="BLAST did not create an output file.",
                )
            )

    except subprocess.CalledProcessError as exc:
        blast_finished_at = datetime.now(timezone.utc)
        elapsed_seconds = time.perf_counter() - blast_started_perf
        logging.error(
            "[%s] BLAST FAILED for locus=%s assembly=%s "
            "(elapsed=%.3f seconds; started=%s)",
            blast_finished_at.isoformat(),
            locus_name,
            assembly_name,
            elapsed_seconds,
            blast_started_at.isoformat(),
        )
        raise RuntimeError(
            f"BLAST failed for locus {locus_name} on {assembly_name}: {exc}"
        ) from exc
    except Exception:
        blast_finished_at = datetime.now(timezone.utc)
        elapsed_seconds = time.perf_counter() - blast_started_perf
        logging.error(
            "[%s] BLAST OUTPUT PARSING FAILED for locus=%s assembly=%s "
            "(elapsed=%.3f seconds; started=%s)",
            blast_finished_at.isoformat(),
            locus_name,
            assembly_name,
            elapsed_seconds,
            blast_started_at.isoformat(),
        )
        raise
    else:
        blast_finished_at = datetime.now(timezone.utc)
        elapsed_seconds = time.perf_counter() - blast_started_perf
        logging.info(
            "[%s] Finished BLAST and parsing for locus=%s assembly=%s "
            "(elapsed=%.3f seconds; started=%s)",
            blast_finished_at.isoformat(),
            locus_name,
            assembly_name,
            elapsed_seconds,
            blast_started_at.isoformat(),
        )
    finally:
        if os.path.exists(blast_output):
            try:
                os.remove(blast_output)
            except OSError as exc:
                logging.warning(
                    "Could not remove temporary BLAST output %s: %s",
                    blast_output,
                    exc,
                )

    return alleles, decisions


def check_for_lock(
    lock_file: Path,
    wait_sec: int = 60,
    stale_after_sec: int | None = 4 * 3600,
) -> None:
    """
    Wait while a schema-writer lock exists.

    This preserves the utility pipeline's protection against reading locus
    FASTAs while another process is updating them. A stale lock can be removed
    after stale_after_sec; active locks are checked at wait_sec intervals.
    """
    lock_file = Path(lock_file)

    while True:
        try:
            stat_result = lock_file.stat()
        except FileNotFoundError:
            return

        if stale_after_sec is not None:
            age = time.time() - stat_result.st_mtime
            if age > stale_after_sec:
                try:
                    lock_file.unlink()
                    print(
                        f"[WARNING] Removed stale schema lock {lock_file} "
                        f"(age {age:.0f} seconds).",
                        file=sys.stderr,
                    )
                    return
                except FileNotFoundError:
                    return

        time.sleep(wait_sec)


def index_loci_fasta(loci_list: list, schema_dir: Path) -> None:
    """
    Compatibility no-op retained for older callers.

    Scheme locus FASTAs are intentionally parsed with Bio.SeqIO and are never
    indexed with pyfaidx. This prevents concurrent samples from creating or
    reading shared, partially written .fai files in the schema directory.
    """
    _ = loci_list, schema_dir
    return


def _emit_pipeline_log(log: object, message: str) -> None:
    """Write a summary message to common pipeline log objects when supplied."""
    if log is None:
        return

    try:
        if hasattr(log, "info") and callable(log.info):
            log.info(message)
            return
        if hasattr(log, "write") and callable(log.write):
            log.write(message + "\n")
            if hasattr(log, "flush") and callable(log.flush):
                log.flush()
    except Exception:
        # Logging must never turn a successful locus call into a failed job.
        print(
            f"[WARNING] Could not write pipeline log message: {message}",
            file=sys.stderr,
        )


def process_single_assembly(
    assembly_path: Path,
    schema_dir: Path,
    output_file: Path,
    log: object,
    max_workers: int,
    # optional knobs (kept optional to avoid changing callers)
    max_stop_extend: int = 0,
    min_cov_ratio: float = 0.70,
    min_identity: float = 90.0,
    region_overlap: float = 0.80,
    stats_file: Path | None = None,
) -> None:
    """
    Run per-locus BLAST for one assembly and write one CDS FASTA.

    This utility-pipeline wrapper preserves the original positional API while
    using the production candidate logic:
    - one representative BLAST hit per physical region;
    - all distinct copies retained;
    - all-allele exact/start/containment checks;
    - no generic frame padding;
    - conservative, opt-in stop extension;
    - start-codon-only normalization to the scheme allele;
    - partial known-allele fragments excluded from the CDS FASTA;
    - a sidecar candidate-decision TSV;
    - the utility pipeline's original temporary-directory and path behavior.
    """
    # Preserve the original pipeline path semantics. Do not resolve symlinks
    # or reinterpret paths relative to a different working directory.
    assembly_path = Path(assembly_path)
    schema_dir = Path(schema_dir)
    output_file = Path(output_file)

    if not assembly_path.exists():
        raise FileNotFoundError(f"Assembly file not found: {assembly_path}")
    if not schema_dir.is_dir():
        raise NotADirectoryError(f"Schema directory not found: {schema_dir}")

    if not 0.0 <= min_cov_ratio <= 1.0:
        raise ValueError("min_cov_ratio must be between 0 and 1")
    if not 0.0 <= min_identity <= 100.0:
        raise ValueError("min_identity must be between 0 and 100")
    if not 0.0 < region_overlap <= 1.0:
        raise ValueError("region_overlap must be greater than 0 and at most 1")
    if max_stop_extend < 0 or max_stop_extend % 3 != 0:
        raise ValueError("max_stop_extend must be a non-negative multiple of 3")

    # Keep the original schema-update lock convention, but do not create or
    # modify any shared pyfaidx indexes in the schema directory.
    check_for_lock(schema_dir / "temp_check.lock")

    loci = [f for f in os.listdir(schema_dir) if f.endswith(".fasta")]
    if not loci:
        raise ValueError(f"No .fasta locus files found in scheme directory: {schema_dir}")

    output_file.parent.mkdir(parents=True, exist_ok=True)
    if stats_file is None:
        stats_file = output_file.with_suffix(".candidate_stats.tsv")
    else:
        stats_file = Path(stats_file)
    stats_file.parent.mkdir(parents=True, exist_ok=True)

    assembly_name = assembly_path.stem
    _emit_pipeline_log(
        log,
        f"Starting BLAST locus calling for {assembly_name} with "
        f"{len(loci)} loci and {max_workers} workers.",
    )
    with tempfile.TemporaryDirectory(
        prefix=f"{fasta_token(assembly_name)}_{os.getpid()}_",
        dir=LOCAL_CWD,
    ) as assembly_tmpdir:
        assembly_index = os.path.join(assembly_tmpdir, "assembly.fai")
        schema_local_dir = Path(os.path.join(assembly_tmpdir, schema_dir.stem))

        try:
            copytree(schema_dir, schema_local_dir, ignore=ignore_patterns('loci_modes'))
        except Exception as exc:
            raise RuntimeError(
                f"Failed to copy schema to local disk: {schema_local_dir}: {exc}"
            ) from exc
  
        try:
            genome = Fasta(str(assembly_path), indexname=assembly_index)
        except Exception as exc:
            raise RuntimeError(
                f"Failed to open assembly FASTA with pyfaidx: {assembly_path}: {exc}"
            ) from exc
    
        try:
            all_alleles = []
            all_decisions = []
    
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = {
                    executor.submit(
                        run_blast_locus,
                        assembly_path,
                        schema_local_dir / locus,
                        genome,
                        assembly_name,
                        min_cov_ratio,
                        min_identity,
                        region_overlap,
                    ): locus
                    for locus in loci
                }
    
                for future in as_completed(futures):
                    locus = futures[future]
                    try:
                        locus_alleles, locus_decisions = future.result()
                    except Exception as exc:
                        # By default: fatal (matches your new script)
                        raise RuntimeError(
                            f"Failed processing locus {locus} for assembly {assembly_path.name}: {exc}"
                        ) from exc
    
                    all_alleles.extend(locus_alleles)
                    all_decisions.extend(locus_decisions)
    
            extracted, extraction_decisions = extract_subsequences(
                genome=genome,
                alleles=all_alleles,
                max_stop_extend=max_stop_extend,
                assembly_name=assembly_name,
            )
            all_decisions.extend(extraction_decisions)
    
            # Write only included candidates to the requested pipeline FASTA.
            with output_file.open("w", encoding="utf-8") as out_fh:
                for header, sequence in extracted:
                    out_fh.write(f"{header}\n{sequence}\n")
    
            write_decision_stats(
                str(stats_file),
                assembly_name,
                all_decisions,
            )
        finally:
            genome.close()

    summary = (
        f"Wrote {len(extracted)} CDS candidates to {output_file}; "
        f"wrote {len(all_decisions)} candidate decisions to {stats_file}."
    )
    print(summary, file=sys.stderr)
    _emit_pipeline_log(log, summary)


