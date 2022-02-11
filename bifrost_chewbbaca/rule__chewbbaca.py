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
        genome = input.genome
        print(genome)
        schemes = params.chewbbaca_schemes
        training_files = params.chewbbaca_training_files # we will remove this, since by default it uses the one in scheme folder
        #schemes = "/local_files/Listeria_monocytogenes_Pasteur_cgMLST_2021-05-31T15"
        #training_files = "/local_files/Listeria_monocytogenes.trn"
        # to do:
        # check if these paths exist, if not make our own scheme and load it
        # need to make sure these folders are mounted properly etc.
        output_dir = output.chewbbaca_results
        os.mkdir(output_dir)
        # copy the contigs file to the output folder, chewbbaca uses a folder containing fastas as input
        shutil.copy(genome, os.path.join(output_dir, sample_name + ".fasta"))
        if species_detection != None:
            species = species_detection["summary"].get("species", None) # currently, provided species will take priority over detected species
        else:
            species = None # in case the category doesnt exist
        #if species not in component["options"]["chewbbaca_current_species"]:
            #species = "Other"
        cmd = f"chewBBACA.py AlleleCall -i {output_dir} -g {schemes} -o {output_dir} --cpu 4"
        print(cmd)
        run_cmd(cmd, log)


    except Exception:
        with open(log.err_file, "w+") as fh:
            fh.write(traceback.format_exc())


rule__chewbbaca(
    snakemake.input,
    snakemake.output,
    snakemake.params,
    snakemake.log)