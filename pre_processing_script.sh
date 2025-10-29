#!/usr/bin/env bash

################################################################################
# BLAST Search and Mapping Script for Sequences
################################################################################
# Description:
#   This pipeline BLASTs input nucleotide sequences against the NCBI nt
#   database (remote) and optionally against SRA runs using blastn_vdb.
#   It then filters hits by identity and coverage and attempts to map
#   sample coordinates for nt accessions via NCBI E-utilities (efetch).
#
# Requirements:
#   - NCBI BLAST+ tools
#   - SRA Toolkit (including blastn_vdb and vdb-validate) [optional]
#   - Internet connection for remote BLAST queries and efetch
################################################################################

# Exit on error, undefined variables, and pipe failures
set -euo pipefail
IFS=$'\n\t'

################################################################################
# Configuration Variables
################################################################################

# Determine script directory and anchor paths
SCRIPT_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

# Input directory (all fasta-like files inside are considered queries)
INPUT_DIR="${SCRIPT_DIR}/01_Input_Sequences"

# Optional: file listing SRA accessions (one per line) to query with blastn_vdb
SRA_ACCESSIONS_FILE="${INPUT_DIR}/sra_accessions.txt"

# BLAST parameters
MIN_IDENTITY=98          # Minimum percent identity threshold
MAX_IDENTITY=100         # Maximum percent identity threshold
: "${MIN_COVERAGE:=0.98}"  # Coverage threshold (fraction of query covered)

# Output directories (numbered by operation order)
OUTPUT_BASE="${SCRIPT_DIR}"
DIR_01_INPUT="${OUTPUT_BASE}/01_Input_Sequences"
DIR_02_NT_BLAST="${OUTPUT_BASE}/02_NT_Database_Results"
DIR_03_SRA_BLAST="${OUTPUT_BASE}/03_SRA_Database_Results"
DIR_04_FILTERED="${OUTPUT_BASE}/04_Filtered_Hits"
DIR_05_COORDINATES="${OUTPUT_BASE}/05_Sample_Coordinates"
TMP_DIR="${OUTPUT_BASE}/.tmp"

################################################################################
# Setup and Validation
################################################################################

echo "Starting BLAST-VDB analysis pipeline..."
echo "========================================"

# Create output directories
mkdir -p "${DIR_01_INPUT}" "${DIR_02_NT_BLAST}" "${DIR_03_SRA_BLAST}" \
         "${DIR_04_FILTERED}" "${DIR_05_COORDINATES}" "${TMP_DIR}"

# Validate input directory exists
if [[ ! -d "${INPUT_DIR}" ]]; then
    echo "Error: Input directory not found: ${INPUT_DIR}"
    exit 1
fi

# Gather query files
shopt -s nullglob
QUERY_FILES=( \
    "${INPUT_DIR}"/*.fa \
    "${INPUT_DIR}"/*.fasta \
    "${INPUT_DIR}"/*.fna \
    "${INPUT_DIR}"/*.fa.gz \
    "${INPUT_DIR}"/*.fasta.gz \
    "${INPUT_DIR}"/*.FA \
    "${INPUT_DIR}"/*.FASTA \
    "${INPUT_DIR}"/*.FNA \
    "${INPUT_DIR}"/*.FA.GZ \
    "${INPUT_DIR}"/*.FASTA.GZ \
)
shopt -u nullglob

if (( ${#QUERY_FILES[@]} == 0 )); then
    echo "Error: No query FASTA files found in ${INPUT_DIR}."
    exit 1
fi

# Snapshot input sequences to organized output directory
for q in "${QUERY_FILES[@]}"; do
    cp -f "$q" "${DIR_01_INPUT}/" 2>/dev/null || true
done

# Check for required tools
command -v blastn >/dev/null 2>&1 || { echo "Error: blastn not found. Install NCBI BLAST+"; exit 1; }
if [[ -s "${SRA_ACCESSIONS_FILE}" ]]; then
    command -v blastn_vdb >/dev/null 2>&1 || echo "Warning: blastn_vdb not found. SRA step will be skipped."
    command -v vdb-validate >/dev/null 2>&1 || echo "Warning: vdb-validate not found. Validation skipped."
fi

################################################################################
# Helper Functions
################################################################################

# Retry wrapper for flaky remote/SRA calls
run_with_retries() {
    local max=${BLAST_MAX_RETRIES:-5}
    local base_delay=${BLAST_BACKOFF_BASE:-3}
    local attempt=1
    while true; do
        if "$@"; then
            return 0
        fi
        local rc=$?
        if (( attempt >= max )); then
            return $rc
        fi
        local sleep_time=$(( base_delay * (2 ** (attempt - 1)) ))
        local jitter=$(( RANDOM % 3 ))
        echo "    Attempt ${attempt} failed (rc=${rc}). Retrying in $((sleep_time + jitter))s..."
        sleep $((sleep_time + jitter))
        attempt=$(( attempt + 1 ))
    done
}

# Decompress .gz files if needed
materialize_query() {
    local src="$1"
    if [[ "$src" =~ \.gz$ ]]; then
        local base="$(basename "$src")"
        local name_no_gz="${base%.gz}"
        local stem="${name_no_gz%%.*}"
        local out="${TMP_DIR}/${stem}.fasta"
        if [[ ! -s "$out" ]]; then
            echo "    Decompressing ${base} -> $(basename "$out")"
            gunzip -c "$src" > "$out"
        fi
        echo "$out"
    else
        echo "$src"
    fi
}

################################################################################
# Step 1: Query NCBI Non-Redundant Nucleotide (nt) Database
################################################################################

echo ""
echo "Step 1: Querying NCBI nt database with megablast (remote)..."
echo "------------------------------------------------------------"

for q in "${QUERY_FILES[@]}"; do
    qbase=$(basename "$q")
    qname_no_gz="${qbase%.gz}"
    qstem="${qname_no_gz%%.*}"
    out_file="${DIR_02_NT_BLAST}/${qstem}_nt_blast.tsv"

    if [[ -s "${out_file}" ]]; then
        echo "  Skipping BLAST nt for ${qbase}: output exists."
        continue
    fi

    q_in="$(materialize_query "$q")"
    echo "  BLAST nt: ${qbase} -> $(basename "${out_file}")"
    run_with_retries blastn \
        -task megablast \
        -query "$q_in" \
        -db nt \
        -remote \
        -outfmt "6 qseqid sacc sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle" \
        -out "$out_file"
done

echo "nt database searches complete."

################################################################################
# Step 2: Query SRA Database (Optional)
################################################################################



################################################################################
# Step 3: Filter Hits by Identity and Alignment Criteria
################################################################################

echo ""
echo "Step 3: Filtering BLAST hits..."
echo "-------------------------------"

filter_blast_hits() {
    local input_file=$1
    local output_file=$2
    awk -v min_id="${MIN_IDENTITY}" -v max_id="${MAX_IDENTITY}" -v min_cov="${MIN_COVERAGE}" '
        {
            qcov = ($6 > 0) ? ($5 / $6) : 0
            if ($4 >= min_id && $4 <= max_id && qcov >= min_cov) {
                print $0
            }
        }
    ' "${input_file}" > "${output_file}"
    echo "  Filtered $(wc -l < "${output_file}") hits from $(basename "${input_file}")"
}

shopt -s nullglob
for f in "${DIR_02_NT_BLAST}"/*.tsv "${DIR_03_SRA_BLAST}"/*.tsv; do
    [[ -e "$f" ]] || continue
    base=$(basename "$f")
    out="${DIR_04_FILTERED}/${base%.tsv}_filtered.tsv"
    if [[ -s "$out" ]]; then
        echo "  Skipping filter for ${base}: output exists."
        continue
    fi
    filter_blast_hits "$f" "$out"
done
shopt -u nullglob

################################################################################
# Step 4: Extract and Map Sample Coordinates
################################################################################
#: << 'OFF'
echo ""
echo "Step 4: Extracting sample coordinates (nt and SRA accessions)..."
echo "----------------------------------------------------------------"

shopt -s nullglob
nt_filtered=("${DIR_04_FILTERED}"/*nt_blast_filtered.tsv)
if (( ${#nt_filtered[@]} > 0 )); then
    cat "${nt_filtered[@]}" | cut -f2 | awk 'length>0' | \
        awk -F'|' '{
            # Extract accession from gi|number|db|accession| format
            if (NF >= 4) {
                # Remove trailing | and version if present
                acc = $4
                gsub(/\|$/, "", acc)
                print acc
            } else {
                print $0
            }
        }' | sort -u > "${DIR_05_COORDINATES}/accessions_nt.txt"
else
    : > "${DIR_05_COORDINATES}/accessions_nt.txt"
fi

# Include SRA accessions if provided
if [[ -s "${SRA_ACCESSIONS_FILE}" ]]; then
    cat "${SRA_ACCESSIONS_FILE}" > "${DIR_05_COORDINATES}/accessions_sra.txt"
else
    : > "${DIR_05_COORDINATES}/accessions_sra.txt"
fi

# Merge both accession lists
if [[ -s "${DIR_05_COORDINATES}/accessions_nt.txt" || -s "${DIR_05_COORDINATES}/accessions_sra.txt" ]]; then
    cat "${DIR_05_COORDINATES}/accessions_nt.txt" "${DIR_05_COORDINATES}/accessions_sra.txt" 2>/dev/null | sort -u > "${DIR_05_COORDINATES}/accessions_all.txt"
else
    : > "${DIR_05_COORDINATES}/accessions_all.txt"
fi
shopt -u nullglob


#OFF
################################################################################
# Step 5: Summary and Cleanup
################################################################################


echo ""
echo "========================================"
echo "Analysis Complete!"
echo "========================================"
echo ""
echo "Results organized in numbered directories:"
echo "  01. Input sequences:      ${DIR_01_INPUT}/"
echo "  02. NT BLAST results:     ${DIR_02_NT_BLAST}/"
echo "  03. SRA BLAST results:    ${DIR_03_SRA_BLAST}/"
echo "  04. Filtered hits:        ${DIR_04_FILTERED}/"
echo "  05. Sample coordinates:   ${DIR_05_COORDINATES}/"
echo ""
if [[ -s "${DIR_05_COORDINATES}/accessions_all.txt" ]]; then
    echo "Total unique accessions processed: $(wc -l < "${DIR_05_COORDINATES}/accessions_all.txt")"
else
    echo "Total unique accessions processed: 0"
fi
if [[ -s "${DIR_05_COORDINATES}/sample_coordinates.tsv" ]]; then
    echo "Samples with coordinates: $(tail -n +2 "${DIR_05_COORDINATES}/sample_coordinates.tsv" | wc -l)"
else
    echo "Samples with coordinates: 0"
fi
echo ""
echo "Pipeline finished successfully."

exit 0
