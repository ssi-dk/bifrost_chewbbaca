import subprocess
import traceback
from bifrostlib import common
from bifrostlib.datahandling import Component, Sample
from bifrostlib.datahandling import SampleComponentReference
from bifrostlib.datahandling import SampleComponent
from bifrostlib.datahandling import Category
from typing import Dict
import os

def run_cmd(command, log):
    with open(log.out_file, "a+") as out, open(log.err_file, "a+") as err:
        command_log_out, command_log_err = subprocess.Popen(command, shell=True).communicate()
        if command_log_err == None:
            command_log_err = ""
        if command_log_out == None:
            command_log_out = ""
        out.write(command_log_out)
        err.write(command_log_err)

def rule__run_cge_resfinder(input: object, output: object, params: object, log: object) -> None:
    try:
        samplecomponent_ref_json = params.samplecomponent_ref_json
        samplecomponent_ref = SampleComponentReference(value=samplecomponent_ref_json)
        samplecomponent = SampleComponent.load(samplecomponent_ref)
        sample = Sample.load(samplecomponent.sample)
        component = Component.load(samplecomponent.component)

        # Variables being used
        species_detection = sample.get_category("species_detection")
        reads = input.reads
        resfinder_db = params.resfinder_db
        pointfinder_db = params.pointfinder_db
        output_dir = output.resfinder_results
        if species_detection != None:
            species = species_detection["summary"].get("species", None) # currently, provided species will take priority over detected species
        else:
            species = None # in case the category doesnt exist
        if species not in component["options"]["resfinder_current_species"]:
            species = "Other"
        species = "\""+species+"\"" # species string must have quotation marks when input to a shell cmd
        if species == '\"Other\"': # this string will be viable input in next resfinder update
            cmd = f"run_resfinder.py -db_res {resfinder_db} -acq -k kma/kma -ifq {reads[0]} {reads[1]} -o {output_dir}"
        else:
            cmd = f"run_resfinder.py -db_res {resfinder_db} -db_point {pointfinder_db} -acq -c -k kma/kma -ifq {reads[0]} {reads[1]} -o {output_dir} --species {species}"
        #print(cmd)
        run_cmd(cmd, log)


    except Exception:
        with open(log.err_file, "w+") as fh:
            fh.write(traceback.format_exc())


rule__run_cge_resfinder(
    snakemake.input,
    snakemake.output,
    snakemake.params,
    snakemake.log)