#- Templated section: start ------------------------------------------------------------------------
import os
import sys
import traceback

from pathlib import Path
from bifrostlib import common
from bifrostlib.datahandling import SampleReference
from bifrostlib.datahandling import Sample
from bifrostlib.datahandling import ComponentReference
from bifrostlib.datahandling import Component
from bifrostlib.datahandling import SampleComponentReference
from bifrostlib.datahandling import SampleComponent
from bifrostlib.datahandling import BioDBReference
from bifrostlib.datahandling import BioDB
from snakemake.io import directory
import datetime

os.umask(0o2)

try:
    sample_ref = SampleReference(_id=config.get('sample_id', None), name=config.get('sample_name', None))
    sample: Sample = Sample.load(sample_ref)
    if sample is None:
        raise Exception("invalid sample passed")

    component_ref = ComponentReference(name=config['component_name'])
    component: Component = Component.load(reference=component_ref)
    if component is None:
        raise Exception("invalid component passed")

    samplecomponent_ref = SampleComponentReference(
        name=SampleComponentReference.name_generator(sample.to_reference(), component.to_reference())
    )
    samplecomponent = SampleComponent.load(samplecomponent_ref)
    if samplecomponent is None:
        samplecomponent = SampleComponent(
            sample_reference=sample.to_reference(),
            component_reference=component.to_reference()
        )

    species_detection = sample.get_category("species_detection")
    species = species_detection["summary"].get("species", None)
    species_sp = species.split()[0]
    print(f"Species is {species} and species_sp is {species_sp}")

    common.set_status_and_save(sample, samplecomponent, "Running")
    
except Exception as error:
    print(traceback.format_exc(), file=sys.stderr)
    raise Exception("failed to set sample, component and/or samplecomponent")

onerror:
    if not samplecomponent.has_requirements():
        common.set_status_and_save(sample, samplecomponent, "Requirements not met")
    if samplecomponent['status'] == "Running":
        common.set_status_and_save(sample, samplecomponent, "Failure")

envvars:
    "BIFROST_INSTALL_DIR",
    "BIFROST_CG_MLST_DIR",
    "CONDA_PREFIX"

rule all:
    input:
        # file is defined by datadump function
        f"{component['name']}/datadump_complete"
    run:
        common.set_status_and_save(sample, samplecomponent, "Success")

rule set_time_start:
    output:
        start_file = f"{component['name']}/time_start.txt"
    run:
        import time
        with open(output.start_file, "w") as fh:
            fh.write(str(time.time()))

rule setup:
    input:
        rules.set_time_start.output.start_file
    output:
        init_file = touch(f"{component['name']}/initialized")
    run:
        samplecomponent['path'] = os.path.join(os.getcwd(), component['name'])
        samplecomponent.save()

rule_name = "check_requirements"
rule check_requirements:
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark"
    input:
        folder = rules.setup.output.init_file,
    output:
        check_file = touch(f"{component['name']}/requirements_met"),
    run:
        if samplecomponent.has_requirements():
            pass
	    
#- Templated section: end --------------------------------------------------------------------------

#* Dynamic section: start **************************************************************************

rule_name = "blast_gene_call"
rule blast_gene_call:
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark"
    input:
        rules.check_requirements.output.check_file,
        genome = f"{sample['categories']['contigs']['summary']['data']}"
    params:
        samplecomponent_ref_json = samplecomponent.to_reference().json,
        chewbbaca_blastdb = f"{os.environ['BIFROST_CG_MLST_DIR']}/blastdb/",
        chewbbaca_schemes = f"{os.environ['BIFROST_CG_MLST_DIR']}/schemes/",
	chunk_output_dir = f"{component['name']}/blast_gene_call_results/fasta_chunks/",
	log_output_dir = f"{component['name']}/blast_gene_call_results/log/",
	chunk_size = 50,
	num_threads = 6
    output:
        gene_call_results = directory(f"{component['name']}/blast_gene_call_results"),
        gene_calls = f"{component['name']}/blast_gene_call_results/gene_calls.fa",
        gene_call_done = f"{component['name']}/blast_gene_call_done"
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "rule__blast_genecall.py")

rule_name = "run_chewbbaca_on_genome"
rule run_chewbbaca_on_genome:
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark"
    input:
        rules.check_requirements.output.check_file,
        rules.blast_gene_call.output.gene_call_done,
        #genome = f"{sample['categories']['contigs']['summary']['data']}"
        genome = rules.blast_gene_call.output.gene_calls
    output:
        chewbbaca_results = directory(f"{component['name']}/chewbbaca_results"),
        #results_tsv = f"{component['name']}/chewbbaca_results/output/results_alleles.tsv",
        chewbbaca_done = f"{component['name']}/chewbbaca_done"
    params:
        samplecomponent_ref_json = samplecomponent.to_reference().json,
        chewbbaca_schemes = f"{os.environ['BIFROST_CG_MLST_DIR']}/schemes/",
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "rule__chewbbaca.py")

#* Dynamic section: end ****************************************************************************

rule set_time_end:
    input:
        rules.run_chewbbaca_on_genome.output.chewbbaca_done
    output:
        end_file = f"{component['name']}/time_end.txt"
    run:
        import time
        with open(output.end_file, "w") as fh:
            fh.write(str(time.time()))

rule_name = "git_version"
rule git_version:
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark"
    input:
        rules.setup.output.init_file
    output:
        git_hash = f"{component['name']}/git_hash.txt"
    run:
        import subprocess, os

        snake_dir = os.path.dirname(workflow.snakefile)

        # Best effort: get commit hash; if not a git repo, write "-"
        try:
            git_hash = subprocess.check_output(
                ["git", "-C", snake_dir, "rev-parse", "HEAD"],
                stderr=subprocess.STDOUT,
                text=True
            ).strip()
        except Exception as e:
            git_hash = "-"
            os.makedirs(os.path.dirname(log.err_file), exist_ok=True)
            with open(log.err_file, "a") as fh:
                fh.write(f"[git_version] Could not determine git hash from {snake_dir}: {e}\n")

        with open(output.git_hash, "w") as fh:
            fh.write(str(git_hash))

rule dump_info:
    input:
        start_file = rules.set_time_start.output.start_file,
        end_file = rules.set_time_end.output.end_file,
        git_hash = rules.git_version.output.git_hash
    output:
        runtime_flag = touch(f"{component['name']}/runtime_set")
    run:
        import time
        from bifrostlib.datahandling import SampleComponent

        with open(input.start_file) as fh:
            t_start = float(fh.read().strip())
        with open(input.end_file) as fh:
            t_end = float(fh.read().strip())
        with open(input.git_hash) as fh:
            git_hash = str(fh.read().strip())
	
        runtime_minutes = (t_end - t_start) / 60.0
        print(f"runtime in minutes {runtime_minutes}")

        sc = SampleComponent.load(samplecomponent.to_reference())
        sc["time_start"] = datetime.datetime.fromtimestamp(t_start).strftime("%Y-%m-%d %H:%M:%S")
        sc["time_end"] = datetime.datetime.fromtimestamp(t_end).strftime("%Y-%m-%d %H:%M:%S")
        sc["time_running"] = round(runtime_minutes, 3)
        sc["git_hash"] = git_hash
	
        sc.save()


#- Templated section: start ------------------------------------------------------------------------
rule_name = "datadump"
rule datadump:
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark"
    input:
        rules.run_chewbbaca_on_genome.output.chewbbaca_results,
        rules.dump_info.output.runtime_flag
    output:
        complete = rules.all.input
    params:
        samplecomponent_ref_json = samplecomponent.to_reference().json
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "datadump.py")
#- Templated section: end --------------------------------------------------------------------------
