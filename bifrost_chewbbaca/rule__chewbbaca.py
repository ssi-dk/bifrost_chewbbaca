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
        command_log_out, command_log_err = subprocess.Popen(command, shell=True).communicate()
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
        sample_name = sample['name']
        # Variables being used
        species_detection = sample.get_category("species_detection")
        detected_species = species_detection['summary']['detected_species']
        genome = input.genome
        schemes = params.chewbbaca_schemes
        species_scheme_folder_matches = [i for i in os.listdir(schemes) if re.match(detected_species.replace(" ", "_") + ".*", i)]
        species_scheme_folder = [i for i in species_scheme_folder_matches if os.path.isdir(os.path.join(schemes, i))][0]
        species_scheme_path = os.path.join(schemes, species_scheme_folder)
        print(species_scheme_path)
        # copy the contigs file to the output folder, chewbbaca uses a folder containing fastas as input
        output_dir = output.chewbbaca_results
        os.mkdir(output_dir)
        #print(output_dir, species_scheme_folder, genome)
        #output_dir = f"{component['name']}"
        shutil.copy(genome, os.path.join(output_dir, sample_name + ".fasta"))
        print(os.listdir(output_dir))
        print(output.chewbbaca_done)
        cmd = f"yes no | chewBBACA.py AlleleCall -i {output_dir} -g {species_scheme_path} -o {output_dir} --cpu 4"
        print(cmd)
        run_cmd(cmd, log)
        with open(output.chewbbaca_done, "w") as fh:
                fh.write("")


    except Exception:
        with open(log.err_file, "w+") as fh:
            fh.write(traceback.format_exc())


rule__chewbbaca(
    snakemake.input,
    snakemake.output,
    snakemake.params,
    snakemake.log)