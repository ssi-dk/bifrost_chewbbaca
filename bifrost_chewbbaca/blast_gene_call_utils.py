#!/usr/bin/env python
import os
import subprocess
import argparse
import sys
import logging
import pathlib
from datetime import datetime
from pathlib import Path 
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import Counter
from Bio.Seq import Seq
from pyfaidx import Fasta
import multiprocessing
import psutil

# Define valid bases and codon sets (uppercase)
VALID_BASES = set(('A', 'T', 'C', 'G'))
START_CODONS = {"ATG", "TTG", "GTG", "CTG", "ATA", "ATT"}
STOP_CODONS  = {"TAA", "TAG", "TGA"}

# Configure logging to file (INFO and above) and preserve warnings/info
logging.basicConfig(
    filename='extraction.log',
    level=logging.INFO,
    format='%(asctime)s %(levelname)s: %(message)s'
)

# Global variables to track start time and previous log time
start_time = datetime.now()
previous_log_time = None

def log_memory_usage(label="Memory usage",log_file="memory_log.txt"):
    global start_time, previous_log_time  # Use global variables to track time
        
    process = psutil.Process(os.getpid())
    memory_info = process.memory_info()
    memory_in_mb = memory_info.rss / (1024 * 1024)  # Convert bytes to MB
    virtual_in_mb = memory_info.vms / (1024 * 1024)
    current_time = datetime.now()
    #current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # Current timestamp
    
    elapsed_time = current_time - start_time
    elapsed_seconds = elapsed_time.total_seconds()
    formatted_elapsed_time = f"{int(elapsed_seconds // 3600):02}:{int((elapsed_seconds % 3600) // 60):02}:{int(elapsed_seconds % 60):02}"
    
    # Calculate time difference from the previous log entry
    if previous_log_time is None:
        time_diff = "00:00:00"
    else:
        delta = current_time - previous_log_time
        delta_seconds = delta.total_seconds()
        time_diff = f"{int(delta_seconds // 3600):02}:{int((delta_seconds % 3600) // 60):02}:{int(delta_seconds % 60):02}"
                                                
    # Update the previous log time
    previous_log_time = current_time

    # Write the log to the specified file
    with open(log_file, "a", encoding="utf-8") as f:
        f.write(f"{label}\tRAM:{memory_in_mb:.2f}MB\tVIRT:{virtual_in_mb:.2f}MB\tWallclock:{current_time.strftime('%Y-%m-%d %H:%M:%S')}\tTotal:{formatted_elapsed_time}\tDifference:{time_diff}\n")

def is_valid_sequence(seq: str) -> bool:
    """
    Check if sequence has only unambiguous DNA bases (case-insensitive).
    """
    return all(base.upper() in VALID_BASES for base in seq)

def is_valid_cds(seq: str) -> bool:
    """
    Returns True if the sequence starts with a valid start codon,
    ends with a valid stop codon, and is a multiple of 3 in length.
    """
    return (
        seq[:3] in START_CODONS and
        seq[-3:] in STOP_CODONS and
        len(seq) % 3 == 0
    )

def has_internal_stop(seq: str) -> bool:
    """
    Returns True if there is any internal in-frame stop codon (excluding first and last codons).
    """
    for i in range(3, len(seq) - 3, 3):
        codon = seq[i:i+3]
        if codon in STOP_CODONS:
            return True
    return False

def reverse_complement(seq: str) -> str:
    """
    Return the reverse complement of a DNA sequence.
    """
    return str(Seq(seq).reverse_complement())

def reverse_complement_check(seq: str) -> bool:
    """
    Decide if sequence should be reverse-complemented based on start codon presence.
    If the forward sequence doesn't start with a known start codon but its reverse
    complement does, return True.
    """
    seq_up = seq.upper()
    if any(seq_up.startswith(c) for c in START_CODONS):
        return False
    rc = reverse_complement(seq)
    rc_up = rc.upper()
    return any(rc_up.startswith(c) for c in START_CODONS)

def orient_and_frame_fix(seq: str, qstart: int, qend: int, genome: Fasta, seq_id: str, header_raw: str):
    """
    Orient (reverse-complement if needed) and ensure seq length % 3 == 0
    by codon-aware padding or trimming. Returns fixed_seq, new_start, new_end, stats.
    """
    stats = Counter()

    # 1) Orientation correction
    if reverse_complement_check(seq):
        seq = reverse_complement(seq)
        stats['reoriented'] += 1

    # Single uppercase pass
    seq = seq.upper()

    # 2) Frame correction
    length = len(seq)
    mod = length % 3
    if mod != 0:
        needed = 3 - mod
        adjusted = False

        contig_len = len(genome[seq_id])

        # a) Try extending 3' end
        if qend + needed <= contig_len:
            cand = genome[seq_id][qstart:qend+needed].seq.upper()
            if (len(cand) % 3 == 0 and
                cand[:3] in START_CODONS and
                cand[-3:] in STOP_CODONS):
                seq = cand
                qend += needed
                stats['padded_stop'] += 1
                adjusted = True

        # b) Try extending 5' end
        if not adjusted and qstart >= needed:
            cand = genome[seq_id][qstart-needed:qend].seq.upper()
            if (len(cand) % 3 == 0 and
                cand[:3] in START_CODONS and
                cand[-3:] in STOP_CODONS):
                seq = cand
                qstart -= needed
                stats['padded_start'] += 1
                adjusted = True

        # c) Fallback: trim
        if not adjusted:
            seq = seq[:length-mod]
            stats['trimmed'] += 1
            msg = (f"{header_raw}: length {length} not multiple of 3; "
                f"trimmed {mod} bases.")
            print(f"[WARNING] {msg}", file=sys.stderr)
            logging.warning(msg)
        else:
            msg = (f"{header_raw}: length {length} not multiple of 3; "
                f"adjusted by {needed} bases.")
            print(f"[INFO]    {msg}", file=sys.stderr)
            logging.info(msg)

    # 3) Extend for nearby start/stop codons
    has_valid_start = seq[:3] in START_CODONS
    has_valid_stop = seq[-3:] in STOP_CODONS

    # Extend for nearby start codons (check positions 1 to 6)
    if not has_valid_start:
        for i in range(1, 7):
            if qstart - i >= 0:
                candidate = genome[seq_id][qstart - i:qend].seq.upper()
                if candidate[:3] in START_CODONS:
                    seq = candidate
                    qstart -= i
                    stats['start_extended'] = i
                    msg = f"{header_raw}: extended upstream by {i} to capture nearby start codon."
                    print(f"[INFO]    {msg}", file=sys.stderr)
                    logging.info(msg)
                    break

    # Extend for nearby stop codons (check up to 6 bp past end)
    if not has_valid_stop:
        for i in range(1, 7):
            if qend + i <= len(genome[seq_id]):
                candidate = genome[seq_id][qstart:qend+i].seq.upper()
                if candidate[-3:] in STOP_CODONS:
                    seq = candidate
                    qend += i
                    stats['stop_extended'] = i
                    msg = f"{header_raw}: extended downstream by {i} to capture nearby stop codon."
                    print(f"[INFO]    {msg}", file=sys.stderr)
                    logging.info(msg)
                    break

    return seq, qstart, qend, stats

def extract_subsequences(genome: Fasta, alleles: list) -> list:
    """
    Iterate over allele regions, perform validation, orientation, framing,
    and return valid (header, seq) pairs.
    """
    extracted = []
    for allele in alleles:
        seq_id, qstart, qend, saccver = allele
        header_raw = f">{seq_id}_{saccver}_{qstart+1}_{qend}"

        # 0) Extract raw subsequence
        subseq = genome[seq_id][qstart:qend].seq

        # 1) Check for ambiguous bases and report details
        invalid_positions = [(i, base) for i, base in enumerate(subseq) if base.upper() not in VALID_BASES]
        if invalid_positions:
            for pos, base in invalid_positions:
                msg = f"{header_raw}: invalid base '{base}' at position {pos+1}"
                print(f"[WARNING] {msg}", file=sys.stderr)
                logging.warning(msg)
            msg = f"{header_raw} contains {len(invalid_positions)} ambiguous or invalid base(s)."
            print(f"[WARNING] {msg}", file=sys.stderr)
            logging.warning(msg)


        # 2) Orient and frame-fix
        fixed, new_start, new_end, stats = orient_and_frame_fix(
            subseq, qstart, qend, genome, seq_id, header_raw)

        # 3) Codon checks
        has_valid_start = any(fixed.startswith(c) for c in START_CODONS)
        has_valid_stop = fixed[-3:] in STOP_CODONS

        # Extended check: nearby start (within first 6 bp)
        nearby_start = any(fixed[i:i+3] in START_CODONS for i in range(1, 7))
        nearby_stop = any(fixed[-i-3:-i] in STOP_CODONS for i in range(1, 7))

        if not has_valid_start:
            if nearby_start:
                msg = f"{header_raw} has no valid start at position 1 but one nearby (within first 6 bp)."
            else:
                msg = f"{header_raw} missing valid start codon."
            print(f"[WARNING] {msg}", file=sys.stderr)
            logging.warning(msg)

        if not has_valid_stop:
            if nearby_stop:
                msg = f"{header_raw} has no valid stop at final 3 bp but one nearby (within last 6 bp)."
            else:
                msg = f"{header_raw} missing valid stop codon."
            print(f"[WARNING] {msg}", file=sys.stderr)
            logging.warning(msg)

        # Check for internal in-frame stop codons
        if has_internal_stop(fixed):
            msg = f"{header_raw}: contains internal in-frame stop codon, likely truncated or pseudogene."
            print(f"[WARNING] {msg}", file=sys.stderr)
            logging.warning(msg)
            continue  # Skip this sequence

        # 4) Final validation and append if valid CDS
        if is_valid_cds(fixed):
            final_hdr = f">{seq_id}_{saccver}_{new_start+1}_{new_end}"
            extracted.append((final_hdr, fixed))
        else:
            msg = f"{header_raw}: skipped due to invalid CDS (start/stop missing or length not divisible by 3)."
            print(f"[WARNING] {msg}", file=sys.stderr)
            logging.warning(msg)
    return extracted

def parse_blast_output(blast_output: str, genome: Fasta) -> list:
    """
    Parse BLAST tabular output and select perfect or best fallback hits.
    Returns alleles as list of tuples: (qaccver, qstart, qend, saccver).
    """
    perfect_hits = []
    fallback = None
    with open(blast_output) as fh:
        for line in fh:
            cols = line.strip().split("\t")
            if len(cols) < 13:
                continue
            rec = {
                'qaccver': cols[0],
                'saccver': cols[1],
                'slen': int(cols[2]),
                'pident': float(cols[3]),
                'length': int(cols[4]),
                'mismatch': int(cols[5]),
                'gapopen': int(cols[6]),
                'qstart': int(cols[7]),
                'qend': int(cols[8]),
                'sstart': int(cols[9]),
                'send': int(cols[10]),
                'evalue': float(cols[11]),
                'bitscore': float(cols[12])
            }
            if rec['pident'] == 100.0 and rec['length'] == rec['slen']:
                perfect_hits.append(rec)
            else:
                if (fallback is None or
                    rec['length'] > fallback['length'] or
                    (rec['length'] == fallback['length'] and rec['pident'] > fallback['pident']) or
                    (rec['length'] == fallback['length'] and rec['pident'] == fallback['pident'] and rec['bitscore'] > fallback['bitscore'])):
                    fallback = rec
    hits = perfect_hits if perfect_hits else ([fallback] if fallback else [])
    alleles = []
    for hit in hits:
        sstart, send = hit['sstart'], hit['send']
        qstart, qend = hit['qstart'], hit['qend']
        slen = hit['slen']
        # Normalize coordinates
        start, end = min(sstart, send), max(sstart, send)
        if start != 1:
            qstart += start - 1
        if end != slen:
            qend += slen - end
            qend = min(qend, len(genome[hit['qaccver']]))
        qstart = max(0, qstart - 1)
        qend = max(qstart, qend)
        alleles.append((hit['qaccver'], qstart, qend, hit['saccver']))
    return alleles

def run_blast_locus(assembly_path: Path, locus_path: Path, genome: Fasta, assembly_name: str) -> list:
    """
    Run BLAST for a single locus against an assembly and return allele regions.
    """
    locus_name = locus_path.name.replace('.fasta', '')
    blast_output = f"blast_{assembly_name}_{locus_name}.txt"
    cmd = [
        'blastn', '-query', assembly_path,
        '-subject', locus_path,
        '-out', blast_output,
        '-outfmt', '6 qaccver saccver slen pident length mismatch gapopen qstart qend sstart send evalue bitscore',
        '-max_target_seqs', '1000'
    ]
    subprocess.run(cmd, check=True)
    alleles = []
    if os.path.exists(blast_output):
        alleles = parse_blast_output(blast_output, genome)
        os.remove(blast_output)
    return alleles

def process_single_assembly(assembly_path: Path, schema_dir: Path, output_file: Path, log: object, max_workers: int):
    if not assembly_path.exists():
        raise FileNotFoundError(f"Assembly file not found: {assembly_path}")

    genome = Fasta(assembly_path)

    loci = [f for f in os.listdir(schema_dir) if f.endswith('.fasta')]

    all_alleles = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                run_blast_locus,
                assembly_path,
                schema_dir/locus,
                genome,
                Path(assembly_path.name).stem
            ): locus for locus in loci
        }
        for future in as_completed(futures):
            all_alleles.extend(future.result())

    extracted = extract_subsequences(genome, all_alleles)

    with open(output_file, 'w') as out_f:
        for hdr, seq in extracted:
            out_f.write(f"{hdr}\n{seq}\n")

    print(f"Wrote {len(extracted)} alleles to {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Run per-locus BLAST in parallel and extract valid CDS subsequences."
    )
    parser.add_argument("--scheme", required=True, help="Directory containing locus FASTA files")
    parser.add_argument("--fa", required=True, help="Directory or file of assembly FASTA sequences")
    parser.add_argument("--out", required=True, help="Directory to store the output files")
    parser.add_argument("--max_workers", type=int, default=multiprocessing.cpu_count(), help="Max concurrent BLASTs")
    parser.add_argument("--samples_file", help="Path to a text file listing one SAMPLE ID per line")

    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)


    args.fa = os.path.abspath(args.fa)
    args.out = os.path.abspath(args.out)
    if args.samples_file:
        args.samples_file = os.path.abspath(args.samples_file)

    # Determine assemblies list
    valid_ext = ('.fa', '.fna', '.fas', '.fasta', '.fa.gz', '.fna.gz', '.fasta.gz')

    if os.path.isdir(args.fa):
        assembly_paths = [
            os.path.join(args.fa, f)
            for f in os.listdir(args.fa)
            if f.lower().endswith(valid_ext)
        ]
    elif os.path.isfile(args.fa):
        assembly_paths = [args.fa]
    else:
        raise ValueError(f"{args.fa!r} is not a file or directory")

    for assembly_path in assembly_paths:
        assembly_name = os.path.splitext(os.path.basename(assembly_path))[0]
        print(f"Processing assembly: {assembly_name}")
        process_single_assembly(
            os.path.basename(assembly_path),
            args.scheme,
            args.out,
            os.path.dirname(assembly_path),
            args.max_workers
        )

if __name__ == "__main__":
    main()
