#!/usr/bin/env bash
set -euo pipefail

# This setup script provisions a conda/mamba environment with BLAST+, SRA Toolkit,
# and helper utilities. Run inside WSL for best results on Windows.

ENV_NAME="VDB_ENV"
CHANNELS=("-c" "conda-forge" "-c" "bioconda")
PKGS=(
    python
    biopython
    blast
    entrez-direct
    seqkit
    sra-tools
    pandas
    gawk jq parallel wget curl
)

have_cmd() { command -v "$1" >/dev/null 2>&1; }

if ! have_cmd conda; then
    echo "Error: conda not found. Install Miniconda/Anaconda or micromamba, then re-run." >&2
    exit 1
fi

# Ensure mamba is available for faster solves
if ! have_cmd mamba; then
    echo "Installing mamba into base..."
    conda install -n base -c conda-forge mamba -y
fi

# Create env if missing
if ! conda env list | grep -q "^${ENV_NAME}\b"; then
    echo "Creating environment: ${ENV_NAME}"
    mamba create -y -n "${ENV_NAME}" ${CHANNELS[@]} ${PKGS[@]}
else
    echo "Environment ${ENV_NAME} already exists. Installing/upgrading packages..."
    mamba install -y -n "${ENV_NAME}" ${CHANNELS[@]} ${PKGS[@]}
fi

echo "\nEnvironment ready. Activate with:"
echo "  conda activate ${ENV_NAME}"
echo "\nVerify tools:"
echo "  blastn -version && blastn_vdb -h | head -n 1 && efetch -help | head -n 1"
