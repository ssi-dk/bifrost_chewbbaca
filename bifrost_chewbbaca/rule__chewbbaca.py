import subprocess
import traceback
from bifrostlib.datahandling import Component, Sample
from bifrostlib.datahandling import SampleComponentReference
from bifrostlib.datahandling import SampleComponent
from pathlib import Path


def run_cmd(command, log, input=None):
    with open(log.out_file, "a+", encoding="utf-8") as out, open(
        log.err_file, "a+", encoding="utf-8"
    ) as err:
        process = subprocess.Popen(
            command,
            shell=False,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            encoding="utf-8",
        )
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
        sample_name = sample["name"]
        # Variables being used
        # resources_dir = component['resources']['schemes']
        species_detection = sample.get_category("species_detection")
        detected_species = species_detection["summary"]["detected_species"]
        try:
            schema = Schema(
                species_name=detected_species,
                schema_home=params.chewbbaca_schemes,
                mapping=component["options"]["chewbbaca_species_mapping"]["schema"],
                log=log
            )
        except KeyError:
            run_cmd(f"touch {component['name']}/no_cgmlst_species_DB", log)
        allelecall = ChewbbacaAlleleCall(
            sample_name, schema, input.genome, output.chewbbaca_results, log
        )
        allelecall.run()
        with open(output.chewbbaca_done, "w", encoding="utf-8") as fh:
            fh.write("")
    except Exception:
        with open(log.err_file, "w+", encoding="utf-8") as fh:
            fh.write(traceback.format_exc())
        raise


class ChewbbacaAlleleCall:
    def __init__(self, sample_name, schema, genome, output, log):
        self.schema = schema
        self.sample_name = sample_name
        self.genome = Path(genome)
        self.output = Path(output)
        self.outputdir = self.output / "output"
        self.log = log
        self.prepare_run()

    def prepare_run(self):
        self.inputdir = self.output / "input"
        self.inputdir.mkdir(parents=True)
        Path(self.inputdir / (self.sample_name + ".fasta")).symlink_to(self.genome.absolute())

    def run(self):
        cmd = [
            "chewBBACA.py",
            "AlleleCall",
            "-i", self.inputdir,
            "-g", self.schema.path(),
            "-o", self.outputdir,
            "--cpu", "4",
            "--pm", "meta",
            "--cds"
        ]
        run_cmd(cmd, self.log, input="no\n")


class Schema:
    def __init__(self, species_name, schema_home, mapping, log) -> None:
        self.species_name = species_name
        self.schema_home = Path(schema_home)
        self.mapping = mapping
        self.log = log
        self.schema_dir = None

    def path(self):
        self.schema_dir = Path(
            self.schema_home, f"{self.mapping[self.species_name]}"
        )
        return self.schema_dir



rule__chewbbaca(snakemake.input, snakemake.output, snakemake.params, snakemake.log)
