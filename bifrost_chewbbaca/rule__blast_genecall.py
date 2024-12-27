import os
import subprocess
import traceback
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from bifrostlib.datahandling import Component, Sample
from bifrostlib.datahandling import SampleComponentReference
from bifrostlib.datahandling import SampleComponent
from pathlib import Path

def run_blastn_and_parse(query_fa, db, assembly_sequences,log):
    """
    Runs BLASTN and processes the output on the fly to keep the best hit per locus-contig pair
    based on sorting criteria, while retaining all hits that are 100% identity and cover the
    full length for rare cases of the same locus appearing twice in one contig.
    """
    blastn_cmd = [
        'blastn',
        '-query', query_fa,
        '-db', db,
        '-outfmt', '6 qaccver saccver slen pident length mismatch gapopen qstart qend sstart send evalue bitscore',
        '-num_threads', '6',
        '-subject_besthit',
        '-max_target_seqs', '2000000',
        '-perc_identity', '90',
        '-max_hsps', '5'
    ]
    print(f"BLAST GENE COMMAND {blastn_cmd}")
    # Dictionary to store the best hit per (locus, contig)
    best_hits = {}
    # List to store all full-coverage, 100%-identity hits
    full_coverage_hits = []

    # Run BLASTN and stream output
    blast_error_file = "blast_error_file"
    with subprocess.Popen(blastn_cmd, stdout=subprocess.PIPE, text=True, stderr=open(blast_error_file, "w+")) as proc:
        stdout,err = proc.communicate()
        if proc.returncode != 0:
            raise RuntimeError(f"Command {' '.join([str(x) for x in blastn_cmd])} failed with code {proc.returncode}.\nCheck the logs in {blast_error_file}")
        
        line_no = 0
        
        for line in stdout.splitlines():
        #for line in stdout:
            # Parse the BLAST output line into a dictionary
            #print(f"DEBUG: {line} and {line.strip().split('\t')}")  # Add this to debug
            cols = line.strip().split("\t")
            line_no = line_no + 1 

            if len(cols) != 13:
                print(f"WARNING: Skipping malformed line: {line.strip()}")
                continue
            
            try:
                record = {
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
                    'bitscore': float(cols[12]),
                }
            except (ValueError, IndexError) as e:
                print(f"ERROR: Skipping line due to parsing error: {line.strip()} ({e})")
                continue
            
            print(f"line number is {line_no}")
            # Extract locus from the subject accession (saccver)
            locus = "_".join(record['saccver'].split("_")[:-1])
            key = (locus, record['qaccver'])  # Locus-contig combination

            # On-the-fly sorting: Keep the best hit per locus-contig pair
            if key in best_hits:
                # Compare current hit with the stored best hit
                best_hit = best_hits[key]
                if (
                    record['pident'] > best_hit['pident'] or  # Higher percent identity
                    (record['pident'] == best_hit['pident'] and record['bitscore'] > best_hit['bitscore']) or  # Higher bitscore
                    (record['pident'] == best_hit['pident'] and record['bitscore'] == best_hit['bitscore'] and record['evalue'] < best_hit['evalue'])  # Lower e-value
                ):
                    best_hits[key] = record
            else:
                # Store this hit as the best for the locus-contig combination
                best_hits[key] = record

            # Collect full-coverage, 100%-identity hits
            if record['length'] == record['slen'] and record['pident'] == 100.0:
                full_coverage_hits.append(record)

    # Collect the final list of hits
    final_hits = list(best_hits.values())

    # Include all full-coverage hits not already in the final list
    for hit in full_coverage_hits:
        if hit not in final_hits:
            final_hits.append(hit)

    # Process the final hits into allele format
    alleles = []
    for hit in final_hits:
        if hit['length'] / hit['slen'] >= 0.6:
            alleles.append(process_hit(hit, assembly_sequences))

    return alleles


def rule__blast_genecall(input: object, output: object, params: object, log: object) -> None:
    try:
        samplecomponent_ref_json = params.samplecomponent_ref_json
        samplecomponent_ref = SampleComponentReference(value=samplecomponent_ref_json)
        samplecomponent = SampleComponent.load(samplecomponent_ref)
        sample = Sample.load(samplecomponent.sample)
        component = Component.load(samplecomponent.component)
        sample_name = sample["name"]
        # Variables being used
        # resources_dir = component['resources']['schemes']
        species_detection = sample.get_category("species_detection")
        detected_species = species_detection["summary"]["detected_species"]

        os.makedirs(output.gene_call_results, exist_ok=True)
        # process_single_assembly(
        #     assembly_file=input.genome, 
        #     db=Path(params.chewbbaca_blastdb)/component["options"]["chewbbaca_species_mapping"]['blastdb'][detected_species],
        #     output_dir=output.gene_call_results, combined_dir, assemblies_dir):
        process_single_assembly(
            assembly_path=input.genome, 
            db=Path(params.chewbbaca_blastdb)/component["options"]["chewbbaca_species_mapping"]['blastdb'][detected_species],
            output_file=output.gene_calls,log=log)

        # process_loci_parallel(
        #     component["options"]["chewbbaca_species_mapping"]['blastdb'][detected_species],
        #     input.genome,
        #     output.gene_call_results)
        with open(output.gene_call_done, "w", encoding="utf-8") as fh:
            fh.write("")
    except Exception:
        with open(log.err_file, "w+", encoding="utf-8") as fh:
            fh.write(traceback.format_exc())
        raise


def process_hit(hit, assembly_sequences):
    """
    Adjusts the start and end positions of a hit and ensures it fits within contig limits.
    Returns a tuple with adjusted hit details.
    """
    alleles=[]
    qaccver, saccver, slen, qlen, sstart, send, qstart, qend = (
        hit['qaccver'], hit['saccver'], int(hit['slen']), int(hit['length']),
        int(hit['sstart']), int(hit['send']), int(hit['qstart']), int(hit['qend'])
    )

    # Ensure that the smallest position is taken as start and the largest as end
    subject_start = min(sstart, send)
    subject_end = max(sstart, send)

    # Adjust start and end points to ensure they're within contig limits
    contig_length = len(assembly_sequences[qaccver])

    if subject_start != 1:
        qstart += subject_start - 1  # Adjust qstart if the alignment doesn't start at position 1
    if subject_end != slen:
        qend += slen - subject_end  # Adjust qend if the alignment doesn't end at the full length
        qend = min(qend, contig_length)

    # Convert qstart to 0-based for BED format and adjust qstart/qend if necessary
    qstart -= 1
    qstart = min(qstart, qend)
    qend = max(qstart, qend)

    return qaccver, qstart, qend, saccver, slen, qlen



def reverse_complement(seq):
    """
    Returns the reverse complement of a sequence.
    """
    return str(Seq(seq).reverse_complement())


def reverse_complement_check(seq):
    """
    Checks if the sequence should be reverse complemented based on its starting pattern.
    """
    return seq.lower().startswith(('tta', 'tca', 'cta'))


def extract_subsequences(fasta_sequences, alleles):
    """
    Extracts subsequences from the given FASTA sequences based on the alleles' positions.
    Returns a list of tuples (header, subsequence).
    """
    extracted_sequences = []
    for allele in alleles:
        seq_id, qstart, qend, saccver = allele[0], allele[1], allele[2], allele[3]
        
        if seq_id in fasta_sequences:
            subseq = fasta_sequences[seq_id][qstart:qend]

            # Check if the subsequence needs to be reverse complemented
            if reverse_complement_check(subseq):
                subseq = reverse_complement(subseq)  # Reverse complement if necessary
            
            header = f">{seq_id}_{saccver}_{qstart + 1}_{qend}"
            extracted_sequences.append((header, subseq))
    
    return extracted_sequences


def read_fasta(file_path):
    """
    Reads a FASTA file and returns a dictionary of sequences.
    """
    return {record.id: str(record.seq) for record in SeqIO.parse(file_path, "fasta")}


def process_single_assembly(assembly_path, db, output_file, log):
    """
    Processes a single assembly against the specified database and writes the combined alleles to a single file.
    """
    assembly_name = os.path.basename(assembly_path).replace('.fasta', '').replace('.fa', '')

    # Read the assembly sequences into memory once per assembly
    fasta_sequences = read_fasta(assembly_path)

    # Run BLAST and parse the output
    alleles = run_blastn_and_parse(assembly_path, db, fasta_sequences, log=log)
    # blast_output = os.path.join(output_dir, f'blast_{assembly_name}.out')
    # run_blastn(assembly_path, db, blast_output)
    
    # # Parse and extract alleles based on the BLAST results
    # alleles = parse_blast_output(blast_output, fasta_sequences)
    extracted_sequences = extract_subsequences(fasta_sequences, alleles)
    
    # Write extracted sequences to the combined FASTA file
    with open(output_file, 'w') as combined_file:
        for header, subseq in extracted_sequences:
            combined_file.write(f"{header}\n{subseq}\n")

rule__blast_genecall(snakemake.input, snakemake.output, snakemake.params, snakemake.log)
