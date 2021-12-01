from bifrostlib import common
from bifrostlib.datahandling import Sample
from bifrostlib.datahandling import SampleComponentReference
from bifrostlib.datahandling import SampleComponent
from bifrostlib.datahandling import Category
from typing import Dict
import os
import json
import re

def extract_resistance(resistance: Category, results: Dict, component_name: str) -> None:
    file_name = "resfinder_results/std_format.json"
    file_key = common.json_key_cleaner(file_name)
    file_path = os.path.join(component_name, file_name)
    with open(file_path) as input:
        results_json = json.load(input)
    results[file_key] = results_json
    phenotypes = results_json['phenotypes']
    seq_regions = results_json['seq_regions']
    seq_variations = results_json['seq_variations']
    dbs = results_json['databases']
    resistance['database_versions'] = {dbs[i]['database_name']:dbs[i]['database_version'] for i in dbs.keys()}
    # collect items for resistant phenotypes
    for phenotype_key in phenotypes.keys():
        phenotype = phenotypes[phenotype_key]
        if phenotype['amr_resistant'] == True:
            amr_classes = phenotype['amr_classes']
            seq_region_keys = phenotype['seq_regions']
            seq_variation_keys = phenotype['seq_variations']
            summary_dict = {
                "name":phenotype_key,
                "seq_regions": ",".join([seq_regions[i]['name'] for i in seq_region_keys]),
                "seq_variations": ",".join([seq_variations[i]['name'] for i in seq_variation_keys]),
                "amr_classes" : ",".join(amr_classes)
            }
            resistance['summary']['phenotypes'].append(summary_dict)

            report_dict = {
                "name":phenotype_key,
                "amr_classes":",".join(amr_classes),
                "seq_regions/variations":[]
            }
            for seq_region in seq_region_keys:
                seq_region_dict = {
                    "name":seq_regions[seq_region]['name'],
                    "coverage":seq_regions[seq_region]['coverage'],
                    "identity":seq_regions[seq_region]['identity'],
                    "type":seq_regions[seq_region]['type'],
                    "seq_var":""
                }
                report_dict['seq_regions/variations'].append(seq_region_dict)
            for seq_variation in seq_variation_keys:
                seq_regions_variant = seq_variations[seq_variation]['seq_regions']
                for seq_region in seq_regions_variant:
                    seq_region_dict = {
                        "name":seq_regions[seq_region]['name'],
                        "coverage":seq_regions[seq_region]['coverage'],
                        "identity":seq_regions[seq_region]['identity'],
                        "type":seq_regions[seq_region]['type'],
                        "seq_var":seq_variations[seq_variation]['seq_var']
                    }
                    report_dict['seq_regions/variations'].append(seq_region_dict)
            resistance['report']['data'].append(report_dict)

def seq_variation_gene_object(gene_list, resfinder_json):
    for gene in gene_list:
        key = resfinder_json['genes'][gene]['key']
        name = resfinder_json['genes'][gene]['name']
        coverage = resfinder_json['genes'][gene]['coverage']
        identity = resfinder_json['genes'][gene]['identity']
        return {"key":key, "name":name, "coverage":coverage, "identity":identity}

def datadump(samplecomponent_ref_json: Dict):
    samplecomponent_ref = SampleComponentReference(value=samplecomponent_ref_json)
    samplecomponent = SampleComponent.load(samplecomponent_ref)
    sample = Sample.load(samplecomponent.sample)
    resistance = samplecomponent.get_category("resistance")
    #print(resistance) # it's the appending that's duplicated because resistance is not none
    #if resistance is None:
    resistance = Category(value={
            "name": "resistance",
            "component": {"id": samplecomponent["component"]["_id"], "name": samplecomponent["component"]["name"]},
            "summary": {"phenotypes":[]},
            "report": {"data":[]}
        }
    )
    extract_resistance(resistance, samplecomponent["results"], samplecomponent["component"]["name"])
    samplecomponent.set_category(resistance)
    sample_category = sample.get_category("resistance")
    if sample_category == None:
        sample.set_category(resistance)
    else:
        current_category_version = extract_digits_from_component_version(resistance['component']['name'])
        sample_category_version = extract_digits_from_component_version(sample_category['component']['name'])
        print(current_category_version, sample_category_version)
        if current_category_version >= sample_category_version:
            sample.set_category(resistance)
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
