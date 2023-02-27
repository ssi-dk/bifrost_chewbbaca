from argparse import Namespace
import pytest
from bifrostlib import common
from bifrostlib import datahandling
from bifrostlib import database_interface
from bifrostlib.datahandling import ComponentReference
from bifrostlib.datahandling import Component
from bifrostlib.datahandling import SampleReference
from bifrostlib.datahandling import Sample
from bifrostlib.datahandling import RunReference
from bifrostlib.datahandling import Run
from bifrost_chewbbaca import launcher
import pymongo
import os, sys
import shutil


@pytest.fixture
def test_connection():
    assert datahandling.has_a_database_connection()
    assert "TEST" in os.environ['BIFROST_DB_KEY'].upper()  # A very basic piece of protection ensuring the word test is in the DB

class TestBifrostchewBBACA:
    component_name = "chewbbaca__v1.0.5"
    component_name = component_name + "__" ## Removed placeholder, since it doesn't get propagated to the component.
    current_dir = os.getcwd()
    test_dir = "/bifrost/test_data/output/test__chewbbaca/"
    json_entries = [
        {
            "_id": {"$oid": "000000000000000000000001"}, 
            "name": "SRR2094561", 
            "components": [], 
            "categories": {
                "contigs": {
                    "summary": {
                        "data": "/bifrost/test_data/samples/SAMN03853109.contigs.fa"
                    }
                },
                "species_detection": {
                    "summary": {
                        "detected_species": "Salmonella enterica"
                    }
                }
            }
        }
    ]
    bson_entries = [database_interface.json_to_bson(i) for i in json_entries]

    @classmethod
    def setup_class(cls):
        with pymongo.MongoClient(os.environ['BIFROST_DB_KEY']) as client:
            db = client.get_database()
            cls.clear_all_collections(db)
            col = db["samples"]
            col.insert_many(cls.bson_entries)
            launcher.initialize()
            os.chdir(cls.current_dir)

    @classmethod
    def teardown_class(cls):
        with pymongo.MongoClient(os.environ['BIFROST_DB_KEY']) as client:
            db = client.get_database()
            #cls.clear_all_collections(db)

    @staticmethod
    def clear_all_collections(db):
        db.drop_collection("components")
        db.drop_collection("hosts")
        db.drop_collection("run_components")
        db.drop_collection("runs")
        db.drop_collection("sample_components")
        db.drop_collection("samples")

    def test_info(self):
        launcher.run_pipeline(["--info"])

    def test_help(self):
        launcher.run_pipeline(["--help"])

    def test_pipeline(self):
        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)

        os.mkdir(self.test_dir)
        test_args = [
            "--sample_name", "SRR2094561",
            "--outdir", self.test_dir
        ]
        launcher.main(args=test_args)
        assert os.path.isfile(f"{self.test_dir}/{self.component_name}/datadump_complete")
        shutil.rmtree(self.test_dir)
        assert not os.path.isdir(f"{self.test_dir}/{self.component_name}")

    def test_db_output(self):
        with pymongo.MongoClient(os.environ['BIFROST_DB_KEY']) as client:
            db = client.get_database()
            samples = db['samples']
            sample_data = next(samples.find({}))
            assert len(sample_data['categories']['cgmlst']['report']['loci']) == 2807
            assert len(sample_data['categories']['cgmlst']['report']['alleles']) == 2807
