

ENV_NAME=$1


SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd $SCRIPT_DIR # avoiding small edge case where bashrc sourcing changes your directory

function exit_function() {
  echo "to rerun use the command:"
  echo "bash -i $SCRIPT_DIR/custom_install.sh $ENV_NAME"
  exit 1
}

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh


if ! git clone https://github.com/ssi-dk/chewBBACA.git
then
    echo >&2 "git clone failed"
    exit_function
fi

cd chewBBACA

if ! (conda env list | grep "$ENV_NAME")
then
 echo "conda environment specified is not found"
 exit_function
else
 conda activate $ENV_NAME
fi

if ! python3 -m pip install --user .
then
    echo >&2 "pip install failed"
    exit_function
fi
