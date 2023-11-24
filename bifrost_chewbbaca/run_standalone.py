import os
import shutil
import pathlib
import subprocess
import sys
from dataclasses import dataclass, asdict
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
    sample_dir = f"{bifrost_install_dir}/bifrost/test_data/samples/"

    @staticmethod
    def clear_all_collections(db):
        db.drop_collection("components")
        db.drop_collection("hosts")
        db.drop_collection("run_components")
        db.drop_collection("runs")
        db.drop_collection("sample_components")
        db.drop_collection("samples")

    def run_pipeline(self):
        os.chdir(self.bifrost_install_dir)
        input_dir = pathlib.Path(self.bifrost_install_dir, 'bifrost', 'test_data', 'samples')
        assert(input_dir.exists())
        child: pathlib.Path
        sample_count = 0
        with pymongo.MongoClient(os.environ["BIFROST_DB_KEY"]) as client:
            db = client.get_database()
            self.clear_all_collections(db)
        for child in input_dir.iterdir():
            if child.is_file() and child.name.endswith('.fasta'):
                sample_count += 1
                if os.path.isdir(self.output_dir):
                    shutil.rmtree(self.output_dir)
                os.makedirs(self.output_dir)
                print()
                print(f"Processing fasta file: {child.name}")
                sample_name = child.name[:-6]
                sample = Sample(
                    name=sample_name,
                    components=list(),
                    categories={
                        'contigs': Category(summary={"data": self.sample_dir + child.name}),
                        'species_detection': Category(summary={"detected_species": "Salmonella enterica"})
                    }
                )
                sample_dict = asdict(sample)
                print("Initial sample document:")
                print(sample_dict)
                with pymongo.MongoClient(os.environ["BIFROST_DB_KEY"]) as client:
                    db = client.get_database()
                    col = db["samples"]
                    bson_entry = database_interface.json_to_bson(sample_dict)
                    col.insert_one(bson_entry)

                test_args = ["--sample_name", sample_name, "--outdir", self.output_dir]

                # This is the important line.
                launcher.main(args=test_args)
                assert (
                    os.path.exists(f"{self.output_dir}/{self.component_name}/datadump_complete")
                    == True
                )
        print(f"Processed {str(sample_count)} samples")
        print("Will now run mongoexport")
        command = \
            f"mongoexport  --db bifrost_test_db --collection samples --pretty --out {self.output_dir}/mongoexport.json"
        process: subprocess.Popen = subprocess.Popen(
            command, stdout=sys.stdout, stderr=sys.stderr, shell=True
        )
        process.communicate()


if __name__ == '__main__':
    start_time = datetime.now()
    instance = BifrostchewBBACA()
    instance.run_pipeline()
    run_time = datetime.now() - start_time
    print(f"Run took {run_time.seconds} seconds")