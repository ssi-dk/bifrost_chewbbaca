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
import psutil
from datetime import datetime
import time
from Bio.SeqRecord import SeqRecord

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

#def split_fasta_in_memory(fasta_path,run_mode,chunk_size=50,length_threshold=1000000,splits=5):
def split_fasta_in_memory(fasta_path,chunk_output_dir,chunk_size=50):
    """
    Splits a FASTA file into chunks of roughly `chunk_size` contigs in memory.
    Returns a list of lists, where each inner list contains FASTA entries.
    """
    
    all_records = list(SeqIO.parse(fasta_path, "fasta"))

    # If the number of records is less than or equal to chunk size, return a single chunk
    if len(all_records) <= chunk_size:
        print(f"the query_fa records are not seperated into chunks")
        return [all_records]  # Wrap in a list to maintain consistency

    contig_chunks = [all_records[i:i + chunk_size] for i in range(0, len(all_records), chunk_size)]

    contig_log_file = os.path.join(chunk_output_dir,"log.txt")
    with open(contig_log_file, "w", encoding="utf-8") as f:
        for idx, chunk in enumerate(contig_chunks):
            f.write(f"Chunk {idx + 1}/{len(contig_chunks)} contains {len(chunk)} contigs:\n")
            for record in chunk:
                f.write(f"  Contig ID: {record.id}, Length: {len(record.seq)}\n")
    print(f"Contig logfile is present at {contig_log_file}")

    return contig_chunks

def split_contigs(query_fa,chunk_output_dir,run_mode,length_threshold,splits):
    #os.makedirs(chunk_output_dir, exist_ok=True) 
    
    fasta_sequences = {}
    original_fasta_sequences = read_fasta(query_fa)
    
    for contig_id, sequence in original_fasta_sequences.items():
        if run_mode == "assembly" and len(sequence) >= length_threshold:
            segment_size = len(sequence) // splits
            # for each segment generated we need to add those as well before seperating them into chunks
            for segment in range(splits):
                start = segment * segment_size + 1
                if segment < splits - 1:
                    end = (segment + 1) * segment_size
                else:
                    end = len(sequence)
                print(f"start {start} and end {end}")
                new_contig_id = f"{contig_id}_segment_{segment+1}"
                fasta_sequences[new_contig_id] = sequence[start:end]
                print(f"\t seperated original {contig_id} : {len(sequence)} into new {new_contig_id} : {len(sequence[start:end])}")
        else:
            fasta_sequences[contig_id] = sequence
            #print(f"Processed contig {contig_id} with length {len(sequence)}")
    
    records = [SeqRecord(Seq(seq), id=contig_id, description="") for contig_id, seq in fasta_sequences.items()]
    
    output_split_fasta = os.path.join(chunk_output_dir,"fasta_contig_split.fasta")

    with open(output_split_fasta, "w") as fasta_file:
        SeqIO.write(records, fasta_file, "fasta")

    print(f"Saved split FASTA sequences to: {output_split_fasta}")

    return output_split_fasta

def run_blastn_and_parse(query_fa, db,log,chunk_output_dir,log_output_dir,chunk_size,num_threads,run_mode,length_threshold,splits):
    """
    Runs BLASTN and processes the output on the fly to keep the best hit per locus-contig pair
    based on sorting criteria, while retaining all hits that are 100% identity and cover the
    full length for rare cases of the same locus appearing twice in one contig.
    """

    #os.makedirs(log_output_dir, exist_ok=True)
    log_file = os.path.join(log_output_dir,"memory_log.txt")
    #print(f"Memory logfile is present at {log_file}")
    
    #os.makedirs(chunk_output_dir, exist_ok=True)
    #output_split_fasta = split_contigs(query_fa,chunk_output_dir,run_mode,length_threshold,splits):
   
    contig_chunks = split_fasta_in_memory(query_fa,chunk_output_dir,chunk_size)

    # Dictionary to store the best hit per (locus, contig)
    best_hits = {}
    full_coverage_hits = []
  
    # Run BLASTN and stream output
    blast_error_file = "blast_error_file"

    log_memory_usage("Memory before launching BLAST subprocess", log_file=log_file)

    for idx, chunk in enumerate(contig_chunks):
        chunk_size = len(chunk)
        print(f"Processing chunk {idx + 1}/{len(contig_chunks)} with {chunk_size} contigs")
        log_memory_usage(f"Memory before BLAST subprocess for chunk {idx + 1} with {chunk_size} contigs", log_file=log_file)
        
        # Convert the chunk back to a FASTA formatted string
        chunk_fasta = "".join(f">{rec.id}\n{str(rec.seq)}\n" for rec in chunk)

        # Ensure chunk is not empty
        if not chunk_fasta.strip():
            print(f"WARNING: Skipping empty chunk {idx + 1}")
            continue
        
        if len(contig_chunks) > 1:
            output_fasta_file = os.path.join(chunk_output_dir, f"chunk_{idx + 1}.fasta")
            with open(output_fasta_file, "w", encoding="utf-8") as fasta_file:
                fasta_file.write(chunk_fasta)
                print(f"FASTA chunk written to {output_fasta_file}")
        
        blastn_cmd = [
            'blastn',
            '-query', '-',  # Use standard input for query
            '-db', db,
            '-outfmt', '6 qaccver saccver slen pident length mismatch gapopen qstart qend sstart send evalue bitscore',
            '-num_threads', str(num_threads),
            '-subject_besthit',
            '-max_target_seqs', '2000000',
            '-perc_identity', '90',
            '-max_hsps', '5'
        ]

        #print(f"Running BLAST command for chunk {idx + 1}/{len(contig_chunks)} \n {' '.join(map(str, blastn_cmd))}")

        with subprocess.Popen(
                blastn_cmd, 
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=open(blast_error_file, "w+"),
                text=True
        ) as proc:

            print(f"Sending input for chunk {idx + 1} to BLAST")
            stdout, err = proc.communicate(input=chunk_fasta) #stdout,_= proc.communicate(input=chunk_fasta)

            log_memory_usage(f"Memory after BLAST subprocess for chunk {idx + 1} with {chunk_size} contigs", log_file=log_file)

            print(f"Finished BLAST for chunk {idx + 1}, processing results...")
            
            if proc.returncode != 0:                       
                if proc.returncode == -9:
                    raise RuntimeError(f"BLAST subprocess was killed by the system (exit code -9) while processing chunk {idx + 1} ")
                else:
                    raise RuntimeError(f"BLAST subprocess failed for chunk {idx + 1} with code {proc.returncode}.")

            #if proc.returncode != 0:
            #    raise RuntimeError(f"Command {' '.join([str(x) for x in blastn_cmd])} failed for cunk {idx + 1} with code {proc.returncode}.\nCheck the logs in {blast_error_file}")
        
            line_no = 0
        
            for line in stdout.splitlines():
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
            
                if line_no % 100000 == 0:
                    print(f"Processed {line_no} lines of chunk {idx + 1} BLAST output")

                # Extract locus from the subject accession (saccver)
                locus = "_".join(record['saccver'].split("_")[:-1])
                key = (locus, record['qaccver'])  # Locus-contig combination

                # On-the-fly sorting: Keep the best hit per locus-contig pair
                if key in best_hits:
                    # Compare current hit with the stored best hit
                    best_hit = best_hits[key]
                    if (record['pident'] > best_hit['pident'] or
                        (record['pident'] == best_hit['pident'] and record['bitscore'] > best_hit['bitscore']) or
                        (record['pident'] == best_hit['pident'] and record['bitscore'] == best_hit['bitscore'] and record['evalue'] < best_hit['evalue'])):
                        best_hits[key] = record
                else:
                    #Store this hit as the best for the locus-contig combination
                    best_hits[key] = record
                            
                # Collect full-coverage, 100%-identity hits
                if record['length'] == record['slen'] and record['pident'] == 100.0:
                    full_coverage_hits.append(record)

    log_memory_usage("After processing BLAST output", log_file=log_file)

    # Collect the final list of hits
    final_hits = list(best_hits.values())

    # Include all full-coverage hits not already in the final list
    for hit in full_coverage_hits:
        if hit not in final_hits:
            final_hits.append(hit)

    #alleles.append(process_hit(hit, assembly_sequences))

    # Process the final hits into allele format
    alleles = []
    for hit in final_hits:
        if hit['length'] / hit['slen'] >= 0.6:
            alleles.append(process_hit(hit, query_fa))

    log_memory_usage("After processing final hits", log_file=log_file)

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
            output_file=output.gene_calls,
            log=log,
            chunk_output_dir=params.chunk_output_dir,
            log_output_dir=params.log_output_dir,
            chunk_size=params.chunk_size,
            num_threads=params.num_threads,
            run_mode=params.run_mode,
            length_threshold=params.length_threshold,
            splits=params.splits
        )

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


def process_hit(hit, fasta_sequences_path):
    """
    Adjusts the start and end positions of a hit and ensures it fits within contig limits.
    Returns a tuple with adjusted hit details.
    """
    
    fasta_sequences = read_fasta(fasta_sequences_path)
    
    alleles=[]
    qaccver, saccver, slen, qlen, sstart, send, qstart, qend = (
        hit['qaccver'], hit['saccver'], int(hit['slen']), int(hit['length']),
        int(hit['sstart']), int(hit['send']), int(hit['qstart']), int(hit['qend'])
    )

    # Ensure that the smallest position is taken as start and the largest as end
    subject_start = min(sstart, send)
    subject_end = max(sstart, send)

    # Adjust start and end points to ensure they're within contig limits
    contig_length = len(fasta_sequences[qaccver])

    if subject_start != 1:
        qstart += subject_start - 1  # Adjust qstart if the alignment doesn't start at position 1
    if subject_end != slen:
        qend += slen - subject_end  # Adjust qend if the alignment doesn't end at the full length
        qend = min(qend, contig_length)

    # Convert qstart to 0-based for BED format and adjust qstart/qend if necessary
    qstart -= 1
    query_start = min(qstart, qend)
    query_end = max(qstart, qend)

    return qaccver, query_start, query_end, saccver, slen, qlen


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


def extract_subsequences(fasta_sequences_path, alleles):
    """
    Extracts subsequences from the given FASTA sequences based on the alleles' positions.
    Returns a list of tuples (header, subsequence).
    """
    fasta_sequences = read_fasta(fasta_sequences_path)


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

def process_single_assembly(assembly_path, db, output_file, log,chunk_output_dir,log_output_dir,chunk_size, num_threads,run_mode,length_threshold,splits):
    """
    Processes a single assembly against the specified database and writes the combined alleles to a single file.
    """
    #assembly_name = os.path.basename(assembly_path).replace('.fasta', '').replace('.fa', '')

    os.makedirs(log_output_dir, exist_ok=True)

    os.makedirs(chunk_output_dir, exist_ok=True)
    output_split_fasta = split_contigs(assembly_path,chunk_output_dir,run_mode,length_threshold,splits)

    # Read the assembly sequences into memory once per assembly
    #original_fasta_sequences = read_fasta(assembly_path) #its in a dictionary
 
    # Run BLAST and parse the output
    alleles = run_blastn_and_parse(
        query_fa=output_split_fasta,
        db=db,
        log=log,
        chunk_output_dir=chunk_output_dir,
        log_output_dir=log_output_dir,
        chunk_size=chunk_size,
        num_threads=num_threads,
        run_mode=run_mode,
        length_threshold=length_threshold,
        splits=splits
    )
    
    extracted_sequences = extract_subsequences(output_split_fasta, alleles)
    
    # Write extracted sequences to the combined FASTA file
    with open(output_file, 'w') as combined_file:
        for header, subseq in extracted_sequences:
            combined_file.write(f"{header}\n{subseq}\n")
    
rule__blast_genecall(snakemake.input, snakemake.output, snakemake.params, snakemake.log)
