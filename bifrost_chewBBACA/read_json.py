import os
import json
file_name = "cge_resfinder__v2_2_3__d98c13b/resfinder_results/std_format_under_development.json"
#file_key = common.json_key_cleaner(file_name)
#file_path = os.path.join(component_name, file_name)
with open(file_name) as input:
    results_json = json.load(input)
#results["asdf"] = results_json
phenotypes = results_json['phenotypes']
genes = results_json['genes']
# collect phenotypes for genes
phen_to_gene_map = {phen:[] for phen in phenotypes.keys()}
for gene in genes.keys():
    gene_phenotype = genes[gene]["phenotypes"]
    if len(gene_phenotype) > 0:
        for phen in gene_phenotype:
                phen_to_gene_map[phen].append(gene)
print(phen_to_gene_map)
for phenotype in phenotypes.keys():
    if phenotypes[phenotype]['resistant']==True:
        phen_dict = {"phenotype":phenotype}
        phen_dict['seq_variations'] = phenotypes[phenotype]['seq_variations']
        phen_dict['genes'] = phen_to_gene_map[phenotype]
        print(phen_dict)