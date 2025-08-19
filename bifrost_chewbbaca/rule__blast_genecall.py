import os
import traceback
from Bio import SeqIO
from Bio.Seq import Seq
#import pandas as pd
from bifrostlib.datahandling import Component, Sample
from bifrostlib.datahandling import SampleComponentReference
from bifrostlib.datahandling import SampleComponent
from pathlib import Path
from blast_gene_call_utils import (
                                   process_single_assembly,
                                   )
#import psutil
#from datetime import datetime
#import time


# def split_fasta_in_memory(fasta_path, chunk_size=50):
#     """
#     Splits a FASTA file into chunks of roughly `chunk_size` contigs in memory.
#     Returns a list of lists, where each inner list contains FASTA entries.
#     """
    
#     all_records = list(SeqIO.parse(fasta_path, "fasta"))

#     # If the number of records is less than or equal to chunk size, return a single chunk
#     if len(all_records) <= chunk_size:
#         print(f"the query_fa records are not seperated into chunks")
#         return [all_records]  # Wrap in a list to maintain consistency

#     chunks = [all_records[i:i + chunk_size] for i in range(0, len(all_records), chunk_size)]
                                               
#     return chunks

# def run_blastn_and_parse(query_fa, db, assembly_sequences,log,chunk_output_dir,log_output_dir,chunk_size, num_threads):
#     """
#     Runs BLASTN and processes the output on the fly to keep the best hit per locus-contig pair
#     based on sorting criteria, while retaining all hits that are 100% identity and cover the
#     full length for rare cases of the same locus appearing twice in one contig.
#     """

#     os.makedirs(log_output_dir, exist_ok=True)
#     log_file = os.path.join(log_output_dir,"memory_log.txt")
#     #print(f"Memory logfile is present at {log_file}")
    
#     os.makedirs(chunk_output_dir, exist_ok=True)
#     contig_chunks = split_fasta_in_memory(query_fa,chunk_size=chunk_size)
#     #print(f"Created fasta chunks are present at {log_file}")
#     contig_log_file = os.path.join(chunk_output_dir,"log.txt")
#     with open(contig_log_file, "w", encoding="utf-8") as f:
#         for idx, chunk in enumerate(contig_chunks):
#             f.write(f"Chunk {idx + 1}/{len(contig_chunks)} contains {len(chunk)} contigs:\n")
#             for record in chunk:
#                 f.write(f"  Contig ID: {record.id}, Length: {len(record.seq)}\n")
#     #print(f"Contig logfile is present at {contig_log_file}")


#     # Dictionary to store the best hit per (locus, contig)
#     best_hits = {}
#     full_coverage_hits = []
  
#     # Run BLASTN and stream output
#     blast_error_file = "blast_error_file"

#     log_memory_usage("Memory before launching BLAST subprocess", log_file=log_file)

#     for idx, chunk in enumerate(contig_chunks):
#         chunk_size = len(chunk)
#         print(f"Processing chunk {idx + 1}/{len(contig_chunks)} with {chunk_size} contigs")
#         log_memory_usage(f"Memory before BLAST subprocess for chunk {idx + 1} with {chunk_size} contigs", log_file=log_file)
        
#         # Convert the chunk back to a FASTA formatted string
#         chunk_fasta = "".join(f">{rec.id}\n{str(rec.seq)}\n" for rec in chunk)

#         # Ensure chunk is not empty
#         if not chunk_fasta.strip():
#             print(f"WARNING: Skipping empty chunk {idx + 1}")
#             continue
        
#         if len(contig_chunks) > 1:
#             output_fasta_file = os.path.join(chunk_output_dir, f"chunk_{idx + 1}.fasta")
#             with open(output_fasta_file, "w", encoding="utf-8") as fasta_file:
#                 fasta_file.write(chunk_fasta)
#                 print(f"FASTA chunk written to {output_fasta_file}")
        
#         blastn_cmd = [
#             'blastn',
#             '-query', '-',  # Use standard input for query
#             '-db', db,
#             '-outfmt', '6 qaccver saccver slen pident length mismatch gapopen qstart qend sstart send evalue bitscore',
#             '-num_threads', str(num_threads),
#             '-subject_besthit',
#             '-max_target_seqs', '2000000',
#             '-perc_identity', '90',
#             '-max_hsps', '5'
#         ]

#         #print(f"Running BLAST command for chunk {idx + 1}/{len(contig_chunks)} \n {' '.join(map(str, blastn_cmd))}")

#         with subprocess.Popen(
#                 blastn_cmd, 
#                 stdin=subprocess.PIPE,
#                 stdout=subprocess.PIPE,
#                 stderr=open(blast_error_file, "w+"),
#                 text=True
#         ) as proc:

#             print(f"Sending input for chunk {idx + 1} to BLAST")
#             stdout, err = proc.communicate(input=chunk_fasta) #stdout,_= proc.communicate(input=chunk_fasta)

#             log_memory_usage(f"Memory after BLAST subprocess for chunk {idx + 1} with {chunk_size} contigs", log_file=log_file)

#             print(f"Finished BLAST for chunk {idx + 1}, processing results...")
            
#             if proc.returncode != 0:
#                 print(f"ERROR: BLAST failed for chunk {idx + 1}")
#                 with open(blast_error_file, "w") as err_file:
#                     err_file.write(stderr)
#                 raise RuntimeError(f"Command {' '.join([str(x) for x in blastn_cmd])} failed for cunk {idx + 1} with code {proc.returncode}.\nCheck the logs in {blast_error_file}")
        
#             line_no = 0
        
#             for line in stdout.splitlines():
#                 cols = line.strip().split("\t")
#                 line_no = line_no + 1 

#                 if len(cols) != 13:
#                     print(f"WARNING: Skipping malformed line: {line.strip()}")
#                     continue
            
#                 try:
#                     record = {
#                         'qaccver': cols[0],
#                         'saccver': cols[1],
#                         'slen': int(cols[2]),
#                         'pident': float(cols[3]),
#                         'length': int(cols[4]),
#                         'mismatch': int(cols[5]),
#                         'gapopen': int(cols[6]),
#                         'qstart': int(cols[7]),
#                         'qend': int(cols[8]),
#                         'sstart': int(cols[9]),
#                         'send': int(cols[10]),
#                         'evalue': float(cols[11]),
#                         'bitscore': float(cols[12]),
#                     }
#                 except (ValueError, IndexError) as e:
#                     print(f"ERROR: Skipping line due to parsing error: {line.strip()} ({e})")
#                     continue
            
#                 #if line_no % 10000 == 0:
#                 #    print(f"Processed {line_no} lines of chunk {idx + 1} BLAST output")

#                 # Extract locus from the subject accession (saccver)
#                 locus = "_".join(record['saccver'].split("_")[:-1])
#                 key = (locus, record['qaccver'])  # Locus-contig combination

#                 # On-the-fly sorting: Keep the best hit per locus-contig pair
#                 if key in best_hits:
#                     # Compare current hit with the stored best hit
#                     best_hit = best_hits[key]
#                     if (record['pident'] > best_hit['pident'] or
#                         (record['pident'] == best_hit['pident'] and record['bitscore'] > best_hit['bitscore']) or
#                         (record['pident'] == best_hit['pident'] and record['bitscore'] == best_hit['bitscore'] and record['evalue'] < best_hit['evalue'])):
#                         best_hits[key] = record
#                 else:
#                     #Store this hit as the best for the locus-contig combination
#                     best_hits[key] = record
                            
#                 # Collect full-coverage, 100%-identity hits
#                 if record['length'] == record['slen'] and record['pident'] == 100.0:
#                     full_coverage_hits.append(record)

#     log_memory_usage("After processing BLAST output", log_file=log_file)

#     # Collect the final list of hits
#     final_hits = list(best_hits.values())

#     # Include all full-coverage hits not already in the final list
#     for hit in full_coverage_hits:
#         if hit not in final_hits:
#             final_hits.append(hit)

#     # Process the final hits into allele format
#     alleles = []
#     for hit in final_hits:
#         if hit['length'] / hit['slen'] >= 0.6:
#             alleles.append(process_hit(hit, assembly_sequences))

#     log_memory_usage("After processing final hits", log_file=log_file)

#     return alleles

# def process_hit(hit, assembly_sequences):
#     """
#     Adjusts the start and end positions of a hit and ensures it fits within contig limits.
#     Returns a tuple with adjusted hit details.
#     """
#     alleles=[]
#     qaccver, saccver, slen, qlen, sstart, send, qstart, qend = (
#         hit['qaccver'], hit['saccver'], int(hit['slen']), int(hit['length']),
#         int(hit['sstart']), int(hit['send']), int(hit['qstart']), int(hit['qend'])
#     )

#     # Ensure that the smallest position is taken as start and the largest as end
#     subject_start = min(sstart, send)
#     subject_end = max(sstart, send)

#     # Adjust start and end points to ensure they're within contig limits
#     contig_length = len(assembly_sequences[qaccver])

#     if subject_start != 1:
#         qstart += subject_start - 1  # Adjust qstart if the alignment doesn't start at position 1
#     if subject_end != slen:
#         qend += slen - subject_end  # Adjust qend if the alignment doesn't end at the full length
#         qend = min(qend, contig_length)

#     # Convert qstart to 0-based for BED format and adjust qstart/qend if necessary
#     qstart -= 1
#     query_start = min(qstart, qend)
#     query_end = max(qstart, qend)

#     return qaccver, query_start, query_end, saccver, slen, qlen

# def reverse_complement_check(seq):
#     """
#     Checks if the sequence should be reverse complemented based on its starting pattern.
#     """
#     return seq.lower().startswith(('tta', 'tca', 'cta'))

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
        
        """
        print(f"species detection {species_detection}\n")
        print(f"detected species {detected_species}\n\n")
        print(f"testing config file {Path(params.chewbbaca_blastdb)}\n")
        print(f"component 1 {component['options']['chewbbaca_species_mapping']}\n")
        print(f"component params {params.chewbbaca_blastdb}\n\n")
        print(f"check test 2 {component['options']['chewbbaca_species_mapping']['blastdb']}\n\n")
        print(f"used database {component['options']['chewbbaca_species_mapping']['blastdb'][detected_species]}\n\n")
        """

        os.makedirs(output.gene_call_results, exist_ok=True)
        process_single_assembly(
            assembly_path=Path(input.genome),
            schema_dir=Path(params.chewbbaca_schemes)/component["options"]["chewbbaca_species_mapping"]['schema'][detected_species],
            # db=Path(params.chewbbaca_blastdb)/component["options"]["chewbbaca_species_mapping"]['blastdb'][detected_species],
            output_file=Path(output.gene_calls),
            log=log,
            max_workers=params.num_threads
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


rule__blast_genecall(snakemake.input, snakemake.output, snakemake.params, snakemake.log)
