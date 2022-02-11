#- Templated section: start ------------------------------------------------------------------------
import os
import sys
import traceback

from bifrostlib import common
from bifrostlib.datahandling import SampleReference
from bifrostlib.datahandling import Sample
from bifrostlib.datahandling import ComponentReference
from bifrostlib.datahandling import Component
from bifrostlib.datahandling import SampleComponentReference
from bifrostlib.datahandling import SampleComponent
os.umask(0o2)

try:
    sample_ref = SampleReference(_id=config.get('sample_id', None), name=config.get('sample_name', None))
    sample:Sample = Sample.load(sample_ref) # schema 2.1
    if sample is None:
        raise Exception("invalid sample passed")
    component_ref = ComponentReference(name=config['component_name'])
    component:Component = Component.load(reference=component_ref) # schema 2.1
    if component is None:
        raise Exception("invalid component passed") # component needs to be in database
    samplecomponent_ref = SampleComponentReference(name=SampleComponentReference.name_generator(sample.to_reference(), component.to_reference()))
    samplecomponent = SampleComponent.load(samplecomponent_ref)
    if samplecomponent is None:
        samplecomponent:SampleComponent = SampleComponent(sample_reference=sample.to_reference(), component_reference=component.to_reference()) # schema 2.1
    common.set_status_and_save(sample, samplecomponent, "Running")

    assemblatron_samplecomponent_field = [i for i in sample['components'] if i['name'].startswith('assemblatron')][0]
    assemblatron_reference = ComponentReference(name=assemblatron_samplecomponent_field['name'])
    assemblatron_samplecomponent_ref = SampleComponentReference(name=SampleComponentReference.name_generator(sample.to_reference(), assemblatron_reference))
    assemblatron_samplecomponent = SampleComponent.load(assemblatron_samplecomponent_ref)
    assemblatron_path = assemblatron_samplecomponent['path']
except Exception as error:
    print(traceback.format_exc(), file=sys.stderr)
    raise Exception("failed to set sample, component and/or samplecomponent")

onerror:
    if not samplecomponent.has_requirements():
        common.set_status_and_save(sample, samplecomponent, "Requirements not met")
    if samplecomponent['status'] == "Running":
        common.set_status_and_save(sample, samplecomponent, "Failure")

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
            #print('asdf')
            with open(output.check_file, "w") as fh:
                fh.write("")

#- Templated section: end --------------------------------------------------------------------------

#* Dynamic section: start **************************************************************************
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
        genome = assemblatron_path
    output:
        chewbbaca_results = directory(f"{component['name']}/chewbbaca_results")
    params:
        samplecomponent_ref_json = samplecomponent.to_reference().json,
        chewbbaca_schemes = component['resources']['schemes'],
        chewbbaca_training_files = component['resources']['training_files']
    #shell:
        #"chewBBACA.py AlleleCall -i /srv/data/FBI/SOFI/EQA/listeria/ -g /srv/data/FBI/SOFI/chewieTest/Listeria_monocytogenes_Pasteur_cgMLST_2021-05-31T15/ --ptf /srv/data/FBI/SOFI/chewieTest/Listeria_monocytogenes.trn -o chewBBACA_output --cpu 5"
        #"run_resfinder.py -db_res {params.resfinder_db} -db_point {params.pointfinder_db} -acq -k kma -ifq {input.reads[0]} {input.reads[1]} -o {output.resfinder_results}"
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
