from argparse import Namespace
import pytest
from bifrostlib import datahandling
from bifrostlib import database_interface
from bifrost_chewbbaca import launcher
import pymongo
import os
import shutil


@pytest.fixture
def test_connection():
    assert datahandling.has_a_database_connection()
    assert (
        "TEST" in os.environ["BIFROST_DB_KEY"].upper()
    )  # A very basic piece of protection ensuring the word test is in the DB


def test_cwd():
    bifrost_install_dir = os.environ["BIFROST_INSTALL_DIR"]
    print(f"bifrost cwd: {bifrost_install_dir}")
    assert bifrost_install_dir != ""


class TestBifrostchewBBACA:
    component_name = "chewbbaca__v1.0.6"

    bifrost_install_dir = os.environ["BIFROST_INSTALL_DIR"]

    test_dir = f"{bifrost_install_dir}/bifrost/test_data/output/test__chewbbaca/"
    r1 = f"{bifrost_install_dir}/bifrost/test_data/samples/SRR2094561.fasta"

    json_entries = [
        {
            "_id": {"$oid": "000000000000000000000001"},
            "name": "SRR2094561",
            "components": [],
            "categories": {
                "contigs": {"summary": {"data": r1}},
                "species_detection": {
                    "summary": {"detected_species": "Salmonella enterica"}
                },
            },
        }
    ]
    bson_entries = [database_interface.json_to_bson(i) for i in json_entries]
    json_biodbs = [
        {
            "_id": {"$oid": "000000000000000000000001"}, 
            "name": "test_cgMLSTschema", 
            "path": "/bifrost/components/bifrost_chewbbaca/resources/test/",
            "summary": {
                "n_alleles": 8,
            }
        }
    ]
    bson_biodbs = [database_interface.json_to_bson(i) for i in json_biodbs]

    @classmethod
    def setup_class(cls):
        with pymongo.MongoClient(os.environ["BIFROST_DB_KEY"]) as client:
            db = client.get_database()
            cls.clear_all_collections(db)
            col = db["samples"]
            col.insert_many(cls.bson_entries)
            col = db["biodbs"]
            col.insert_many(cls.bson_biodbs)
            launcher.initialize()
            os.chdir(cls.bifrost_install_dir)

    @classmethod
    def teardown_class(cls):
        with pymongo.MongoClient(os.environ["BIFROST_DB_KEY"]) as client:
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
        db.drop_collection("biodbs")

    def test_info(self):
        launcher.run_pipeline(["--info"])

    def test_help(self):
        launcher.run_pipeline(["--help"])

    def test_pipeline(self):
        if os.path.isdir(self.test_dir):
            shutil.rmtree(self.test_dir)

        os.makedirs(self.test_dir)
        test_args = ["--sample_name", "SRR2094561", "--outdir", self.test_dir]
        launcher.main(args=test_args)
        assert (
            os.path.exists(f"{self.test_dir}/{self.component_name}/datadump_complete")
            == True
        )
        #shutil.rmtree(self.test_dir)
        assert not os.path.isdir(f"{self.test_dir}/{self.component_name}")

    def test_db_output(self):
        with pymongo.MongoClient(os.environ["BIFROST_DB_KEY"]) as client:
            print(f"databases: {client.list_database_names()}")
            db = client.get_database()
            print(f"collections: {db.list_collection_names()}")
            samples = db["samples"]
            sample_data = next(samples.find({}))
            assert len(sample_data['categories']['cgmlst']['report']['loci']) == 2804
            assert len(sample_data['categories']['cgmlst']['report']['alleles']) == 2804
