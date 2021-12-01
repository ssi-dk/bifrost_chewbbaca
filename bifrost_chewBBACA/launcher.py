#!/usr/bin/env python3
"""
Launcher file for accessing dockerfile commands
"""
import argparse
from io import FileIO, StringIO
import json
import os
import sys
import traceback
import subprocess
from bifrostlib import datahandling
from bifrostlib.datahandling import Component
from bifrostlib.datahandling import ComponentReference
import yaml
import pprint
from typing import List, Dict

from yaml.tokens import StreamEndToken


global COMPONENT


def initialize():
    with open(os.path.join(os.path.dirname(__file__), 'config.yaml')) as fh:
        config: Dict = yaml.load(fh, Loader=yaml.FullLoader)

    if not(datahandling.has_a_database_connection()):
        raise ConnectionError("BIFROST_DB_KEY is not set or other connection error")

    global COMPONENT
    try:
        component_ref = ComponentReference(name=config["name"])
        COMPONENT = Component.load(component_ref)
        if COMPONENT is not None and '_id' in COMPONENT.json:
                return
        else:
            COMPONENT = Component(value=config)
            install_component()

    except Exception as e:
        print(traceback.format_exc(), file=sys.stderr)
    return


def install_component():
    COMPONENT['install']['path'] = os.path.os.getcwd()
    print(f"Installing with path:{COMPONENT['install']['path']}")
    try:
        COMPONENT.save()
        print(f"Done installing")
    except:
        print(traceback.format_exc(), file=sys.stderr)
        sys.exit(0)


class types():
    def file(path):
        if os.path.isfile(path):
            return os.path.abspath(path)
        else:
            raise argparse.ArgumentTypeError(f"{path} #Bad file path")

    def directory(path):
        if os.path.isdir(path):
            return os.path.abspath(path)
        else:
            raise argparse.ArgumentTypeError(f"{path} #Bad directory path")


def parse_and_run(args: List[str]) -> None:
    description: str = (
        f"-Description------------------------------------\n"
        f"{COMPONENT['details']['description']}"
        f"------------------------------------------------\n"
        f"\n"
        f"-Environmental Variables/Defaults---------------\n"
        f"BIFROST_CONFIG_DIR: {os.environ.get('BIFROST_CONFIG_DIR','.')}\n"
        f"BIFROST_RUN_DIR: {os.environ.get('BIFROST_RUN_DIR','.')}\n"
        f"BIFROST_DB_KEY: {os.environ.get('BIFROST_DB_KEY')}\n"
        f"------------------------------------------------\n"
        f"\n"
    )

    # Using two parsers for UX so that install doesn't conflict while all the args still point to the main parser
    basic_parser = argparse.ArgumentParser(add_help=False)
    basic_parser.add_argument(
        '--reinstall',
        action='store_true',
    )
    basic_parser.add_argument(
        '--info',
        action='store_true',
        help='Provides basic information on COMPONENT'
    )

    #Second parser for the arguements related to the program, everything can be set to defaults (or has defaults)
    parser: argparse.ArgumentParser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Show arg values'
    )
    parser.add_argument(
        '-out', '--outdir',
        default=os.environ.get('BIFROST_RUN_DIR', os.getcwd()),
        help='Output directory'
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        '-name', '--sample_name',
        action='store',
        type=str,
        help='Sample name of sample in bifrost, sample has already been added to the bifrost DB'
    )
    group.add_argument(
        '-id', '--sample_id',
        action='store',
        type=str,
        help='Sample ID of sample in bifrost, sample has already been added to the bifrost DB'
    )

    try:
        basic_options, extras = basic_parser.parse_known_args(args)
        if basic_options.reinstall:
            install_component()
            return None
        elif basic_options.info:
            show_info()
            return None
        else:
            pipeline_options, junk = parser.parse_known_args(extras)
            if pipeline_options.debug is True:
                print(pipeline_options)
            run_pipeline(pipeline_options)
    except Exception as e:
        print(traceback.format_exc, file=sys.stderr)


def show_info() -> None:
    pprint.pprint(COMPONENT.json)


def run_pipeline(args: argparse.Namespace) -> None:
    try:
        sample_var = f""
        if args.sample_id is not None:
            sample_var = f"sample_id={args.sample_id}"
        else:
            sample_var = f"sample_name={args.sample_name}"
        command = f"cd {args.outdir}; snakemake --nolock --cores all -s {os.path.join(os.path.dirname(__file__),'pipeline.smk')} --config {sample_var} component_name={COMPONENT['name']}"
        print(command)
        process: subprocess.Popen = subprocess.Popen(
            command,
            stdout=sys.stdout,
            stderr=sys.stderr,
            shell=True
        )
        process.communicate()
    except:
        print(traceback.format_exc())


def main(args = sys.argv):
    initialize()
    parse_and_run(args)


if __name__ == '__main__':
    main()
