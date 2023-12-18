import os
import shutil
import pathlib
import subprocess
import sys
from dataclasses import dataclass, asdict
from datetime import datetime
import random
import string

import pymongo

from bifrostlib import database_interface
from bifrost_chewbbaca import launcher

"""
Loads a new sample into MongoDB and runs ChewBBACA on it.
"""

def rndstr(length):
    return ''.join(random.choices(string.ascii_uppercase + string.digits, k=length))

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
        run_name = rndstr(5)
        with pymongo.MongoClient(os.environ["BIFROST_DB_KEY"]) as client:
            db = client.get_database()
            self.clear_all_collections(db)
        
        # Check that we do not get too long sample names!
        for child in input_dir.iterdir():
            if child.is_file() and (child.name.endswith('.fasta') or child.name.endswith('.fna')):
                sample_name = child.name[:-6]
                if len(sample_name) > 30:
                    print(f"File name too long: {child.name}")
                    sys.exit()

        for child in input_dir.iterdir():
            if child.is_file() and (child.name.endswith('.fasta') or child.name.endswith('.fna')):
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
                        'species_detection': Category(summary={"detected_species": "Salmonella enterica"}),
                        'sample_info':  Category(summary={
                            'institution': 'TEST',
                            'sofi_sequence_id': run_name + '_' + sample_name
                            })
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


if __name__ == '__main__':
    export_dir = pathlib.Path(os.environ["MONGOEXPORT_DIR"])
    assert export_dir.exists()
    start_time = datetime.now()
    instance = BifrostchewBBACA()
    instance.run_pipeline()
    now = datetime.now()
    run_time = now - start_time
    print(f"Run took {run_time.seconds} seconds")
    print("Will now run mongoexport")
    export_filename = now.isoformat(timespec='seconds').replace(':', '') + '.json'
    export_filepath = pathlib.Path(export_dir, export_filename)
    print(f"Export will be saved in {export_filepath}")
    command = \
        f"mongoexport --db bifrost_test_db --collection samples --pretty --out {str(export_filepath)}"
    process: subprocess.Popen = subprocess.Popen(
        command, stdout=sys.stdout, stderr=sys.stderr, shell=True
    )
    process.communicate()