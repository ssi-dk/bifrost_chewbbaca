import subprocess
import traceback
from bifrostlib import common
from bifrostlib.datahandling import Component, Sample
from bifrostlib.datahandling import SampleComponentReference
from bifrostlib.datahandling import SampleComponent
from bifrostlib.datahandling import Category
from typing import Dict
import os
import shutil
import re

def run_cmd(command, log):
    with open(log.out_file, "a+") as out, open(log.err_file, "a+") as err:
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
        command_log_out, command_log_err = process.communicate()
        if command_log_err == None:
            command_log_err = ""
        if command_log_out == None:
            command_log_out = ""
        out.write(command_log_out)
        err.write(command_log_err)


def rule__chewbbaca(input: object, output: object, params: object, log: object) -> None:
    try:
        samplecomponent_ref_json = params.samplecomponent_ref_json
        samplecomponent_ref = SampleComponentReference(value=samplecomponent_ref_json)
        samplecomponent = SampleComponent.load(samplecomponent_ref)
        sample = Sample.load(samplecomponent.sample)
        component = Component.load(samplecomponent.component)
        print(component['options'])
        sample_name = sample['name']
        # Variables being used
        resources_dir = component['resources']['schemes']
        species_detection = sample.get_category("species_detection")
        detected_species = species_detection['summary']['detected_species']
        try:
            species_name, species_id, schema_id = get_species_id(detected_species, component["options"]["chewbbaca_species_mapping"], resources_dir, log)
        except KeyError:
            run_cmd(f"touch {component['name']}/no_cgmlst_species_DB", log)
        genome = input.genome
        schemes = params.chewbbaca_schemes
        print(schemes)
        # copy the contigs file to the output folder, chewbbaca uses a folder containing fastas as input
        output_dir = output.chewbbaca_results
        os.mkdir(output_dir)
        shutil.copy(genome, os.path.join(output_dir, sample_name + ".fasta"))
        cmd = f"yes no | chewBBACA.py AlleleCall -i {output_dir} -g {resources_dir}/{species_name}_{species_id}_{schema_id}/* -o {output_dir} --cpu 4"
        print(cmd)
        run_cmd(cmd, log)
        sync_schema(species_name, species_id, schema_id, resources_dir, log)
        with open(output.chewbbaca_done, "w") as fh:
                fh.write("")


    except Exception:
        with open(log.err_file, "w+") as fh:
            fh.write(traceback.format_exc())

def get_species_id(species_name, mapping, schema_dir, log):
    mapping[species_name]
    fetch_schema(*mapping[species_name], schema_dir, log)
    return mapping[species_name]

def fetch_schema(schema_name, species_id, schema_id, schema_dir, log):
    if not os.path.exists(f"{schema_dir}/{schema_name}_{species_id}_{schema_id}"):
        cmd = f"chewBBACA.py DownloadSchema -sp {species_id} -sc {schema_id} -o {schema_dir}/{schema_name}_{species_id}_{schema_id}"
        run_cmd(cmd,log)
    else:
        sync_schema(schema_name, species_id, schema_id, schema_dir, log)

def sync_schema(schema_name, species_id, schema_id, schema_dir, log):
    cmd = f"chewBBACA.py SyncSchema -sc {schema_dir}/{schema_name}_{species_id}_{schema_id}/*"
    run_cmd(cmd, log)



rule__chewbbaca(
    snakemake.input,
    snakemake.output,
    snakemake.params,
    snakemake.log)