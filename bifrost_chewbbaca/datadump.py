from bifrostlib import common
from bifrostlib.datahandling import Sample
from bifrostlib.datahandling import SampleComponentReference
from bifrostlib.datahandling import SampleComponent
from bifrostlib.datahandling import Category
from typing import Dict
import os
import json
import re
from pathlib import Path

def call_percent(calls):
    return round(len([x for x in calls if x.isdigit() or x.startswith('INF')])/len(calls)*100,2)

def extract_cgmlst(cgmlst: Category, results: Dict, component_name: str) -> None:
    output_folder = Path(component_name, 'chewbbaca_results', 'output')
    # chewbacca output gets thrown into a folder called results_<yearmonthday>someothertext
    #chewbbaca_output_folder = [i for i in os.listdir(output_folder) if re.match("results_[0-9]{6}.*", i)][0]
    file_name = output_folder/ "results_alleles.tsv"
    file_key = common.json_key_cleaner(str(file_name))
    with open(file_name) as input:
        lines = input.readlines()
        lines = [i.strip() for i in lines]
        locus_names = lines[0].split()[1:]
        allele_values = lines[1].split()[1:]
        allele_dict = {locus_names[i]:allele_values[i] for i in range(len(locus_names))}
        cgmlst['summary']['call_percent'] = call_percent(allele_values)
    results[file_key] = allele_dict
    cgmlst['report']['alleles'] = allele_dict
    cgmlst['report']['loci'] = locus_names


def datadump(samplecomponent_ref_json: Dict):
    samplecomponent_ref = SampleComponentReference(value=samplecomponent_ref_json)
    samplecomponent = SampleComponent.load(samplecomponent_ref)
    sample = Sample.load(samplecomponent.sample)
    cgmlst = samplecomponent.get_category("cgmlst")
    if cgmlst is None:
        cgmlst = Category(value={
                "name": "cgmlst",
                "component": {"id": samplecomponent["component"]["_id"], "name": samplecomponent["component"]["name"]},
                "summary": {"call_percent": None},
                "report": {
                    "loci":[],
                    "alleles": {},
                    }
                }
            )
    extract_cgmlst(cgmlst, samplecomponent["results"], samplecomponent["component"]["name"])
    samplecomponent.set_category(cgmlst)
    sample_category = sample.get_category("cgmlst")
    if sample_category == None:
        sample.set_category(cgmlst)
    else:
        current_category_version = extract_digits_from_component_version(cgmlst['component']['name'])
        sample_category_version = extract_digits_from_component_version(sample_category['component']['name'])
        print(current_category_version, sample_category_version)
        if current_category_version >= sample_category_version:
            sample.set_category(cgmlst)
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
