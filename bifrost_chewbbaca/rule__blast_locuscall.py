import os
import traceback
from Bio import SeqIO
from Bio.Seq import Seq
#import pandas as pd
from bifrostlib.datahandling import Component, Sample
from bifrostlib.datahandling import SampleComponentReference
from bifrostlib.datahandling import SampleComponent
from pathlib import Path
from blast_locus_call_utils import (
                                   process_single_assembly,
                                   )

import sys

sys.stdout = open(snakemake.log.out_file, "a")
sys.stderr = open(snakemake.log.err_file, "a")


def rule__blast_locuscall(input: object, output: object, params: object, log: object) -> None:
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
        
        print(f"species detection {species_detection}\n")
        print(f"detected species {detected_species}\n\n")
        print(f"testing config file {Path(params.chewbbaca_blastdb)}\n")
        print(f"component 1 {component['options']['chewbbaca_species_mapping']}\n")
        print(f"component params {params.chewbbaca_blastdb}\n\n")
        print(f"check test 2 {component['options']['chewbbaca_species_mapping']['blastdb']}\n\n")
        print(f"used database {component['options']['chewbbaca_species_mapping']['blastdb'][detected_species]}\n\n")

        os.makedirs(output.locus_call_results, exist_ok=True)
        process_single_assembly(
            assembly_path=Path(input.genome),
            schema_dir=Path(params.chewbbaca_schemes)/component["options"]["chewbbaca_species_mapping"]['schema'][detected_species],
            output_file=Path(output.locus_calls),
            log=log,
            max_workers=6
        )

        with open(output.locus_call_done, "w", encoding="utf-8") as fh:
            fh.write("")
    except Exception:
        with open(log.err_file, "w+", encoding="utf-8") as fh:
            fh.write(traceback.format_exc())
        raise


rule__blast_locuscall(snakemake.input, snakemake.output, snakemake.params, snakemake.log)
