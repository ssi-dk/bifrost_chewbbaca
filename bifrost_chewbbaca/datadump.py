from bifrostlib import common
from bifrostlib.datahandling import Sample
from bifrostlib.datahandling import SampleComponentReference
from bifrostlib.datahandling import SampleComponent
from bifrostlib.datahandling import Category
from typing import Dict
import os
import json
import re

def extract_cgmlst(chewbbaca: Category, results: Dict, component_name: str) -> None:
    output_folder = os.path.join(component_name, 'chewbbaca_results')
    # chewbacca output gets thrown into a folder called results_<yearmonthday>someothertext
    chewbbaca_output_folder = [i for i in os.listdir(output_folder) if re.match("results_[0-9]{6}.*", i)][0]
    file_name = os.path.join("chewbbaca_results", chewbbaca_output_folder, "results_alleles.tsv")
    file_key = common.json_key_cleaner(file_name)
    file_path = os.path.join(component_name, file_name)
    with open(file_path) as input:
        lines = input.readlines()
        lines = [i.strip() for i in lines]
        allele_names = lines[0].split()[1:]
        allele_values = lines[1].split()[1:]
        allele_dict = {allele_names[i]:allele_values[i] for i in range(len(allele_names))}
    results[file_key] = allele_dict
    #chewbbaca['summary']['alleles'] = allele_dict
    chewbbaca['report']['data'].append({"alleles":allele_dict})


def datadump(samplecomponent_ref_json: Dict):
    samplecomponent_ref = SampleComponentReference(value=samplecomponent_ref_json)
    samplecomponent = SampleComponent.load(samplecomponent_ref)
    sample = Sample.load(samplecomponent.sample)
    #chewbbaca = samplecomponent.get_category("chewbbaca")
    #print(resistance) # it's the appending that's duplicated because resistance is not none
    #if resistance is None:
    chewbbaca = Category(value={
            "name": "chewbbaca",
            "component": {"id": samplecomponent["component"]["_id"], "name": samplecomponent["component"]["name"]},
            "summary": {"sequence_type":None},
            "report": {"data":[]}
        }
    )
    extract_cgmlst(chewbbaca, samplecomponent["results"], samplecomponent["component"]["name"])
    samplecomponent.set_category(chewbbaca)
    sample_category = sample.get_category("chewbbaca")
    if sample_category == None:
        sample.set_category(chewbbaca)
    else:
        current_category_version = extract_digits_from_component_version(chewbbaca['component']['name'])
        sample_category_version = extract_digits_from_component_version(sample_category['component']['name'])
        print(current_category_version, sample_category_version)
        if current_category_version >= sample_category_version:
            sample.set_category(chewbbaca)
    common.set_status_and_save(sample, samplecomponent, "Success")
    
    with open(os.path.join(samplecomponent["component"]["name"], "datadump_complete"), "w+") as fh:
        fh.write("done")


def extract_digits_from_component_version(component_str):
    version_re = re.compile(".*__(v.*)__.*")
    version_group = re.match(version_re, component_str).groups()[0]
    version_digits = int("".join([i for i in version_group if i.isdigit()]))
    return version_digits
datadump(
    snakemake.params.samplecomponent_ref_json,
)
