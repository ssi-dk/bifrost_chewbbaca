#!/bin/bash

set -euo pipefail

helpFunction()
{
   echo ""
   echo "Usage:"
   echo "$0 -i LOCAL - for local install"
   echo "$0 -i COMP - for computerome install"
   exit 1 # Exit script after printing help
}

while getopts "i:" opt
do
   case "$opt" in
      i ) parameterI="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$parameterI" ]
then
   echo "-i parameter is empty";
   helpFunction
fi

if [ "$parameterI" != "LOCAL" ] && [ "$parameterI" != "COMP" ]
then
  echo "Wrong argument"
  helpFunction
fi

if [ "$parameterI" == "LOCAL" ]
then
  echo "Starting local install"
fi

if [ "$parameterI" == "COMP" ]
then
  echo "Starting computerome install"
  module load tools computerome_utils/2.0
  module load tools $CONDA_VERSION
  #if $BIFROST_CONDA_PATH is not set then exit with help message
  if [ -z "$BIFROST_CONDA_PATH" ]
  then
    echo "Please set $BIFROST_CONDA_PATH variable to your prefered env install location"
    echo "Example:"
    echo "export BIFROST_CONDA_PATH=/path/to/env/install/location"
    exit 1
  else
    echo -e "\nAdding conda envs_dirs and pkgs_dirs"
    echo "conda config --add envs_dirs $BIFROST_CONDA_PATH/envs"
    echo "conda config --add pkgs_dirs $BIFROST_CONDA_PATH/pkgs"
    conda config --add envs_dirs $BIFROST_CONDA_PATH/envs
    conda config --add pkgs_dirs $BIFROST_CONDA_PATH/pkgs
  fi
fi

# Begin script
if $(conda config --show channels | grep -q "bioconda")
then
  echo "bioconda channel is already added"
else
  conda config --add channels bioconda
  echo "bioconda channel was added, you can remove it after installation with command"
  echo "conda config --remove channels bioconda"
fi

if $(conda config --show channels | grep -q "conda-forge")
then
  echo "conda-forge channel is already added"
else
  conda config --add channels conda-forge
  echo "conda-forge channel was added, you can remove it after installation with command"
  echo "conda config --remove channels conda-forge"
fi

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
REQ_TXT="$SCRIPT_DIR/environment.yml"

CONFIG_YAML_PATH=$(find $SCRIPT_DIR -name "config.yaml")
if test -f "$CONFIG_YAML_PATH";
then
  COMPONENT_NAME=$(grep "display_name:.*." $CONFIG_YAML_PATH | tr " " "\n" | grep -v "display_name:")
  if [ -z "$COMPONENT_NAME" ]
    then
      echo "display_name: in config.yaml should contain component name"
      exit 1
  fi
  COMPONENT_VERSION=$(grep -o "code:.*." $CONFIG_YAML_PATH | tr " " "\n" | grep -v "code:")
  if [ -z "$COMPONENT_VERSION" ]
    then
      echo "code: in config.yaml should contain component version"
      exit 1
  fi
  ENV_NAME=("bifrost_"$COMPONENT_NAME"_"$COMPONENT_VERSION)
else
  echo "Cannot find config.yaml in component folder to form env name"
  exit 1
fi

CUSTOM_INSTALL_PATH=$(find $SCRIPT_DIR -name "custom_install.sh")
#check if env $ENV_NAME already exists
if ! (conda env list | grep "$ENV_NAME")
then
  echo "Environment is not found and will be created"
else
  echo "Environment $ENV_NAME already exists"
  echo -e "\nIf you want to update it, please remove it first"
  echo "conda env remove --name $ENV_NAME"
  if test -f "$CUSTOM_INSTALL_PATH"
  then
    echo -e "\nIf you want to run the custom install part, you can execute:"
    echo "bash -i $CUSTOM_INSTALL_PATH $ENV_NAME"
  fi
  exit 1
fi

#check if environment.yml file exists
if test -f "$REQ_TXT";
then
  echo "Making conda env"
  echo "$ENV_NAME will be created"
  mamba env create --file "$REQ_TXT" --name $ENV_NAME
else
  echo "environment.yml file cannot be found in the script folder"
  exit 1
fi

if ! (conda env list | grep "$ENV_NAME")
then
  if [ "$parameterI" == "COMP" ]
  then
    echo -e "\nRemoving conda envs_dirs and pkgs_dirs"
    echo "conda config --remove envs_dirs $BIFROST_CONDA_PATH/envs"
    echo "conda config --remove pkgs_dirs $BIFROST_CONDA_PATH/pkgs"
    conda config --remove envs_dirs $BIFROST_CONDA_PATH/envs
    conda config --remove pkgs_dirs $BIFROST_CONDA_PATH/pkgs
  fi
  echo "Environment $ENV_NAME was not created"
  echo "Inspect conda error messages"
  exit 1
else
  echo "Environment $ENV_NAME was created"
fi

#check if custom_install.sh file exists and run it
if test -f "$CUSTOM_INSTALL_PATH";
then
  echo -e "\nRunning custom_install.sh"
  bash -i $CUSTOM_INSTALL_PATH $ENV_NAME #-i required for interactive mode to active env
fi
