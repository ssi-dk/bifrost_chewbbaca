import os
import shutil
import pathlib
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime

import pymongo

from bifrostlib import database_interface
from bifrost_chewbbaca import launcher

"""
Loads a new sample into MongoDB and runs ChewBBACA on it.
"""


@dataclass
class Category:
    summary: dict

@dataclass
class Sample:
    name: str
    components: list
    categories: dict

class BifrostchewBBACA:
    component_name = "chewbbaca__v1.0.6"

    bifrost_install_dir = os.environ["BIFROST_INSTALL_DIR"]

    output_dir = f"{bifrost_install_dir}/bifrost/test_data/output/test__chewbbaca/"
    sample_dir = f"{bifrost_install_dir}/bifrost/test_data/samples/SRR2094561.fasta"

    sample_template = {
        "name": "SRR2094561",
        "components": [],
        "categories": {
            "contigs": {"summary": {"data": sample_dir}},
            "species_detection": {
                "summary": {"detected_species": "Salmonella enterica"}
            },
        },
    }

    @staticmethod
    def clear_all_collections(db):
        db.drop_collection("components")
        db.drop_collection("hosts")
        db.drop_collection("run_components")
        db.drop_collection("runs")
        db.drop_collection("sample_components")
        db.drop_collection("samples")

    def test_pipeline(self):
        with pymongo.MongoClient(os.environ["BIFROST_DB_KEY"]) as client:
            db = client.get_database()
            self.clear_all_collections(db)
            col = db["samples"]
            bson_entry = database_interface.json_to_bson(self.sample_template)
            col.insert_one(bson_entry)

        os.chdir(self.bifrost_install_dir)
        if os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)

        os.makedirs(self.output_dir)
        input_dir = pathlib.Path(self.bifrost_install_dir, 'bifrost', 'test_data', 'samples')
        assert(input_dir.exists())
        child: pathlib.Path
        # Currently there can only be one fasta file ind the folder, otherwise test will fail.
        # So the for loop is actually meaningless.
        for child in input_dir.iterdir():
            if child.is_file() and child.name.endswith('.fasta'):
                sample_name = child.name[:-6]
                test_args = ["--sample_name", sample_name, "--outdir", self.output_dir]

                # This is the important line.
                launcher.main(args=test_args)
                assert (
                    os.path.exists(f"{self.output_dir}/{self.component_name}/datadump_complete")
                    == True
                )

                command = \
                    f"mongoexport  --db bifrost_test_db --collection samples --pretty --out {self.output_dir}/mongoexport.json"
                process: subprocess.Popen = subprocess.Popen(
                    command, stdout=sys.stdout, stderr=sys.stderr, shell=True
                )
                process.communicate()


if __name__ == '__main__':
    start_time = datetime.now()
    instance = BifrostchewBBACA()
    instance.test_pipeline()
    run_time = datetime.now() - start_time
    print(f"Run took {run_time.seconds} seconds")