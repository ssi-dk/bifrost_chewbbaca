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

def rule__split_assembly_fasta(input: object, output: object, params: object, log: object) -> None:
    output_records = []

    length_threshold = params.length_threshold
    splits = params.splits
    output_fasta = output.split_genome

    try:
        # Read all sequences in the FASTA file
        for record in SeqIO.parse(input.genome, "fasta"):
            seq_length = len(record.seq)

        if seq_length >= length_threshold:
            # Split the sequence into `params.splits` parts
            chunk_size = seq_length // splits
            print(f"Splitting {record.id} (length: {seq_length}) into {splits} parts of {chunk_size} bases each")

            for i in range(splits):
                start = i * chunk_size
                end = (i + 1) * chunk_size if i < splits - 1 else seq_length  # ensures the final piece of last chunk is also extracted
                chunk_seq = record.seq[start:end]

                # Create new record with modified header
                chunk_record = record[:]
                chunk_record.id = f"{record.id}_part{i+1}"
                chunk_record.description = f"Split from {record.id}, bases {start+1}-{end}"
                chunk_record.seq = chunk_seq

                output_records.append(chunk_record)

        else:
            # Keep the contig as it is
            print(f"Keeping {record.id} (length: {seq_length}) unchanged")
            output_records.append(record)

        # Write all modified sequences to output FASTA
        with open(output_fasta, "w") as f:
            SeqIO.write(output_records, f, "fasta")
    except Exception:
        with open(log.err_file, "w+", encoding="utf-8") as fh:
            fh.write(traceback.format_exc())
        raise

    print(f"\nSplitting completed! Output saved to {output.split_genome}")

rule__split_assembly_fasta(snakemake.input, snakemake.output, snakemake.params, snakemake.log)
