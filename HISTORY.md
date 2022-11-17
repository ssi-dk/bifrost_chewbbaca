# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v0_0_1] - 2022-02-07
Initializing component

## [v1_0_0] - 2022-11-17
Schemas from ChewieNS + tests now working

Added code to rule__chewbbaca.py to fetch cgMLST schemas from central chewieNS server
Fixed warlock version to 1.3.3, since the latest warlock version hangs when validating against bifrost schema.
Changed chewBBACA version to 3.0.0 from github repo.
Disentangled bifrost_chewbbaca component from assemblatron component, and grabs contigs from location indicated in bifrostdb.
Fixed component name in test_simple, so datadump test can pass.
updated miniconda3 to version 4.12.0 in order to be compatible with current version of prodigal
