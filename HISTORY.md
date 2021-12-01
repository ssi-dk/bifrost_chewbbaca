# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v2_2_8] - 2021-02-17
Updates to bifrostlib to fix datetime bug

### Changed
- setup.cfg
## [v2_2_7] - 2021-02-12
### Notes
Update bifrostlib and cleanup dockerfile for potentially not pulling on base
### Changed
- Dockerfile
- setup.py

## [v2_2_6] - 2021-02-12
### Notes
Fix bug missing num of reads in report
### Changed
- bifrost_min_read_check/datadump.py

## [v2_2_5] - 2021-02-11
### Notes
Update dependencies bifrostlib to 2.1.5
### Changed
- setup.py
## [v2_2_2] - 2020-12-17
### Notes
Fixed bump2version to also update test_simple.py
### Changed
- setup.cfg
- /tests/test_simple.py
## [v2_2_2] - 2020-12-17
### Notes
Forgot to update history before bumping, but made sure it worked on github for CI/CD
### Changed
- Dockerfile
- .github/workflows/run_tests.yml
- setup.cfg
  - Better utilizing bump2vesion options
## [v2_2_1] - 2020-12-17
### Notes
Development has adjusted from dev/prod on github to test/prod. Test for testing and prod for live code. Dev is now local and utilizes the whole bifrost folder structure while test/prod only utilizes the submodule structure. Currently tests have not been updated to reflect this change.

### Added
- test_simple.py
- setup.cfg 
  - For bump2version and pytest init

### Changed
Adjusted format of docker to work in context to root folder of project for local development /bifrost/ and to the submodule to prod. New environment test also is created for testing prod ready code on github actions with CI. 
- .dockerignore
- Dockerfile
- docker_build_and_push_to_dockerhub.yml
- test_standard_workflow.yml -> run_tests.yml
Requirements are now for pytest and environmental values, actual required python libraries are included in setup.py and required programs in the Dockerfile
- requirements.txt
- setup.py
Adjust code to reflect new bifrostlib==2.2.0 which has a new data structure and schema for the objects
- /bifrost_min_read_check/\_\_init\_\_.py
- /bifrost_min_read_check/\_\_main\_\_.py
- /bifrost_min_read_check/config.yaml
- /bifrost_min_read_check/launcher.py
- /bifrost_min_read_check/datadump.py
- /bifrost_min_read_check/pipeline.smk
- /bifrost_min_read_check/rule__greater_than_min_reads_check.py
- CHANGELOG.md -> HISTORY.md

### Removed
Development is now local and the dockerfile for compose on the root /bifrost/ folder instead of at each submodule
- docker-compose.dev.yaml
- docker-compose.yaml
Now unnecessary as requirements merged with setup.py
- requirements.dev.txt
- test_1_standard_workflow.py
  - Testing structure adjusted now using test_simple.py

## [2.2.0] - 2020-10-07
### Added
- CHANGELOG.md into repo

### Changed
- Templated files to be common for SSI maintained pipelines. All files use .env and <COMPONENT_NAME>/config.yaml as the source of all information. The config.yaml should be considered the primary source of all information regarding the component and it's settings. The .env file needs to contain the <COMPONENT_NAME> and install specific settings (currently just mongo_db connection). 
  - docker-compose.dev.yaml
  - docker-compose.yaml
  - .env
    - This is being used for both Dockerfile and passing the values into the Docker image, not sure if that has any issues.
  - setup.py
    - This can't use libraries to extract config values, so right now it's hardcoded on what to look for. This can cause some potential issues.
- The following files are also impacted by the changes
  - Dockerfile

### Removed
- Docker-compose files no longer point to a custom env file