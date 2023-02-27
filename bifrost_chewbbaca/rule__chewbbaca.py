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
from pathlib import Path

def run_cmd(command, log, input=None):
    with open(log.out_file, "a+") as out, open(log.err_file, "a+") as err:
        process = subprocess.Popen(command, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
        command_log_out, command_log_err = process.communicate(input=input)
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
            schema = Schema(detected_species,params.chewbbaca_schemes, component["options"]["chewbbaca_species_mapping"],component['resources']['genelists'], log)
        except KeyError:
            run_cmd(f"touch {component['name']}/no_cgmlst_species_DB", log)
        allelecall = ChewbbacaAlleleCall(sample_name, schema, input.genome, output.chewbbaca_results, log)
        allelecall.run()
        schema.sync_schema()
        with open(output.chewbbaca_done, "w") as fh:
                fh.write("")
    except Exception:
        with open(log.err_file, "w+") as fh:
            fh.write(traceback.format_exc())
        raise

class ChewbbacaAlleleCall:
    def __init__(self, sample_name, schema, genome, output, log):
        self.schema = schema
        self.sample_name = sample_name
        self.genome = Path(genome)
        self.output = Path(output)
        self.outputdir = self.output/"output"
        self.log = log
        self.prepare_run()

    def prepare_run(self):
        self.inputdir = self.output/"input"
        self.inputdir.mkdir(parents=True)
        Path(self.inputdir/(self.sample_name+".fasta")).symlink_to(self.genome)

    def run(self):
        cmd = ["chewBBACA.py","AlleleCall","-i",self.inputdir,"-g",self.schema.path(),"-o",self.outputdir,"--cpu", "4"]
        cmd.extend(self.schema.get_genelist())
        run_cmd(cmd, self.log, input="no\n")


class Schema:
    def __init__(self, species_name, schema_home, mapping, genelist_dir, log) -> None:
        self.species_name = species_name
        self.schema_home = Path(schema_home)
        self.mapping = mapping
        self.log = log
        self.schema_dir = None
        self.genelist_dir = Path(genelist_dir)
        self.get_species_id()

    def path(self):
        self.schema_path = Path(self.schema_home, f"{self.schema_name}_{self.species_id}_{self.schema_id}")
        if not self.schema_path.exists():
            self.fetch_schema()
        if self.schema_dir is None:
            self.schema_dir = next(self.schema_path.iterdir())
        return self.schema_dir

    def get_genelist(self):
        ## Returns a genelist option for the chewbbaca command or an empty string if no such list exists.
        genelist = self.genelist_dir / f"{self.species_id}_{self.schema_id}.list"
        if genelist.exists():
            return ["--gl", str(genelist)]
        else:
            return []

    def get_species_id(self):
        self.schema_name, self.species_id, self.schema_id = self.mapping[self.species_name]
        return self.species_id

    def fetch_schema(self):
        cmd = ["chewBBACA.py", "DownloadSchema", "-sp", str(self.species_id), "-sc", str(self.schema_id), "-o" ,self.schema_path]
        run_cmd(cmd,self.log)

    def sync_schema(self):
        cmd = ["chewBBACA.py", "SyncSchema", "-sc", str(self.path())]
        run_cmd(cmd, self.log)



rule__chewbbaca(
    snakemake.input,
    snakemake.output,
    snakemake.params,
    snakemake.log)