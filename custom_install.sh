#!/usr/bin/env bash
# custom_install.sh
# Usage: bash -i custom_install.sh <ENV_NAME>
# Runs inside an interactive shell so "conda activate" works (install.sh already does this).

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <CONDA_ENV_NAME>"
  exit 1
fi

ENV_NAME="$1"

# Require credentials as env vars
: "${GIT_USER:?Please export GIT_USER before running install}"
: "${GIT_TOKEN:?Please export GIT_TOKEN before running install}"

# Activate the env created by install.sh
# (install.sh calls: bash -i custom_install.sh $ENV_NAME, so conda should be available)
conda activate "$ENV_NAME"

# --no-deps to keep dependency resolution within your conda env
pip install --no-deps "git+https://${GIT_USER}:${GIT_TOKEN}@github.com/ssi-dk/chewbbaca_fbi.git"

echo "chewbbaca_fbi installation complete."
