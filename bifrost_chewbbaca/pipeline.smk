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
os.umask(0o2)


try:
    sample_ref = SampleReference(_id=config.get('sample_id', None), name=config.get('sample_name', None))
    sample:Sample = Sample.load(sample_ref) # schema 2.1
    if sample is None:
        raise Exception("invalid sample passed")
    sample_name =sample['name']

    species_detection = sample.get_category("species_detection")
    species = species_detection["summary"].get("species", None)
    species_sp = species.split()[0]
    print(f"Species is {species} and species_sp is {species_sp}")

    sample_info = sample.get_category("sample_info")
    run_mode_sequencing = sample_info["summary"].get("haveReads")
    run_mode_assembly = sample_info["summary"].get("haveAsm")
    mode = None
    print(f"{run_mode_sequencing}")
   
    if run_mode_sequencing == True and run_mode_assembly == False:
        print(f"running sequencing mode {run_mode_sequencing}")
        mode = "sequencing"
        #ruleorder: setup > check_requirements > blast_gene_call > run_chewbbaca_on_genome > datadump > all
    elif run_mode_assembly == True and run_mode_sequencing == False:
        print(f"running assembly mode {run_mode_assembly}")
        mode = "assembly"
	#ruleorder: setup > check_requirements > split_assembly_fasta > blast_gene_call > run_chewbbaca_on_genome > datadump > all
    elif run_mode_sequencing == False and run_mode_assembly == False:
        print(f"failed to run any chewbacca mode")
    elif run_mode_sequencing == True and run_mode_assembly == True:
        mode = "sequencing"
	#ruleorder: setup > check_requirements > blast_gene_call > run_chewbbaca_on_genome > datadump > all
    
    component_ref = ComponentReference(name=config['component_name'])
    component:Component = Component.load(reference=component_ref) # schema 2.1
    if component is None:
        raise Exception("invalid component passed") # component needs to be in database
    samplecomponent_ref = SampleComponentReference(name=SampleComponentReference.name_generator(sample.to_reference(), component.to_reference()))
    samplecomponent = SampleComponent.load(samplecomponent_ref)
    try:
        print(samplecomponent.has_requirements())
    except Exception:
        samplecomponent:SampleComponent = SampleComponent(sample_reference=sample.to_reference(), component_reference=component.to_reference()) # schema 2.1
        print(f"{samplecomponent.json}")
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

rule setup:
    output:
        init_file = touch(temp(f"{component['name']}/initialized")),
    params:
        folder = component['name']
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
        check_file = f"{component['name']}/requirements_met",
    params:
        samplecomponent
    run:
        if samplecomponent.has_requirements():
          if not Path(output.check_file).exists():
            with open(output.check_file, "w") as fh:
                fh.write("")

#rule_name = "split_assembly_fasta"
#rule split_assembly_fasta:
#    message:
#        f"Running step:{rule_name}"
#    log:
#        out_file = f"{component['name']}/log/{rule_name}.out.log",
#        err_file = f"{component['name']}/log/{rule_name}.err.log",
#    benchmark:
#        f"{component['name']}/benchmarks/{rule_name}.benchmark"
#    input:
#        rules.check_requirements.output.check_file,
#        genome = f"{sample['categories']['contigs']['summary']['data']}" 
#    params:
#        log_output_dir = f"{component['name']}/log",
#        length_threshold = 1000000,
#        splits = 5	
#    output:
#        split_genome = f"{component['name']}/split_fasta.fa"
#    script:
#        os.path.join(os.path.dirname(workflow.snakefile), "rule__split_assembly_fasta.py")

#    shell:
#        "echo {output.split_genome} with {params.length_threshold} and {params.splits}" 

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
	chunk_output_dir = f"{component['name']}/blast_gene_call_results/fasta_chunks/",
	log_output_dir = f"{component['name']}/blast_gene_call_results/log/",
	chunk_size = 25,
	num_threads = 6,
	run_mode = mode,
	length_threshold = 1000000,
	splits = 5
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
        #* Dynamic section: start ******************************************************************
        rules.run_chewbbaca_on_genome.output.chewbbaca_results  # Needs to be output of final rule
        #* Dynamic section: end ********************************************************************
    output:
        complete = rules.all.input
    params:
        samplecomponent_ref_json = samplecomponent.to_reference().json
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "datadump.py")
#- Templated section: end --------------------------------------------------------------------------
