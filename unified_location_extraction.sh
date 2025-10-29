#!/usr/bin/env bash

################################################################################
# Unified Location and Metadata Extraction Script
################################################################################
# Description:
#   Comprehensive metadata extraction combining ALL working methods for:
#   - NT (NCBI Nucleotide) database accessions
#   - SRA (Sequence Read Archive) accessions
#
#   Methods integrated:
#   1. GenBank flat file lat_lon field extraction
#   2. BioSample lat_lon via NT accession link
#   3. BioSample lat_lon via SRA accession link
#   4. SRA runinfo for platform/library information
#   5. Comprehensive metadata extraction (country, geo_loc, strain, etc.)
#
#   Output: Single unified TSV file with all available metadata
#
# Requirements:
#   - NCBI EDirect tools (esearch, elink, efetch, xtract)
#   - Internet connection for NCBI queries
#
# Author: Automated Pipeline
# Date: October 29, 2025
################################################################################

set -uo pipefail

################################################################################
# Environment Setup
################################################################################

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate VDB_ENV

# Change to working directory
cd /mnt/c/Users/Mark_Cyril/PIPELINE/Kuya_Elds/BLAST_VDB

################################################################################
# Configuration
################################################################################

OUTPUT_DIR="05_Sample_Coordinates"
FINAL_OUTPUT="${OUTPUT_DIR}/unified_metadata.tsv"
NT_ACCESSIONS="${OUTPUT_DIR}/accessions_nt.txt"
SRA_ACCESSIONS="${OUTPUT_DIR}/accessions_sra.txt"
TMP_DIR="/tmp/blast_vdb_unified_$$"
LOG_FILE="${OUTPUT_DIR}/extraction.log"

# Create directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${TMP_DIR}"

# Cleanup function
cleanup() {
    rm -rf "${TMP_DIR}"
}
trap cleanup EXIT

# Initialize log - redirect stderr to log file, keep stdout for data
exec 2>> "${LOG_FILE}"

################################################################################
# Helper Functions
################################################################################

# Parse lat/lon from various formats and convert to decimal degrees
parse_latlon() {
    local latlon="$1"
    echo "$latlon" | awk '{
        # Handle formats like:
        # "35.99 N 120.42 E"
        # "35.99, 120.42"
        # "-35.99 120.42"
        # "35.99N 120.42E"
        gsub(/[,;]/, " ")
        lat=""; lon=""; lat_sign=1; lon_sign=1
        for(i=1; i<=NF; i++) {
            if ($i ~ /^-?[0-9]+\.?[0-9]*$/) {
                if (lat == "") {
                    lat = $i
                    if (lat < 0) { lat_sign = -1; lat = -lat }
                } else if (lon == "") {
                    lon = $i
                    if (lon < 0) { lon_sign = -1; lon = -lon }
                }
            } else if (tolower($i) ~ /^[ns]$/) {
                if (tolower($i) == "s") lat_sign = -1
            } else if (tolower($i) ~ /^[ew]$/) {
                if (tolower($i) == "w") lon_sign = -1
            }
        }
        if (lat != "" && lon != "") {
            printf "%.6f\t%.6f", lat*lat_sign, lon*lon_sign
        } else {
            print "N/A\tN/A"
        }
    }'
}

################################################################################
# NT Accession Extraction (Comprehensive Multi-Method)
################################################################################

extract_nt_metadata() {
    local acc="$1"
    local tmpfile="${TMP_DIR}/nt_${acc}.gb"
    
    echo "  [NT] Fetching GenBank record..." >&2
    
    # Initialize all variables
    local organism="N/A"
    local biosample="N/A"
    local country="N/A"
    local geo_loc="N/A"
    local latitude="N/A"
    local longitude="N/A"
    local isolate="N/A"
    local strain="N/A"
    local cultivar="N/A"
    local collection_date="N/A"
    local host="N/A"
    local tissue="N/A"
    local platform="N/A"
    local library="N/A"
    local lat_lon=""
    
    # ========================================================================
    # METHOD 1: GenBank Flat File Parsing
    # ========================================================================
    if ! efetch -db nuccore -id "${acc}" -format gb 2>/dev/null > "${tmpfile}"; then
        echo "  ✗ Failed to fetch GenBank record" >&2
        echo -e "${biosample}\t${organism}\t${country}\t${geo_loc}\t${latitude}\t${longitude}\t${isolate}\t${strain}\t${cultivar}\t${collection_date}\t${host}\t${tissue}\t${platform}\t${library}"
        return 1
    fi
    
    # Extract organism
    organism=$(grep -m1 "ORGANISM" "${tmpfile}" | sed 's/.*ORGANISM  *//' | tr -d '\n' | sed 's/[[:space:]]*$//' || echo "N/A")
    
    # Extract all available metadata fields from GenBank
    country=$(grep -i '/country=' "${tmpfile}" | head -n1 | sed -E 's/.*country="?([^"]+)"?.*/\1/' | tr -d '\n' || echo "N/A")
    geo_loc=$(grep -i '/geo_loc_name=' "${tmpfile}" | head -n1 | sed -E 's/.*geo_loc_name="?([^"]+)"?.*/\1/' | tr -d '\n' || echo "N/A")
    lat_lon=$(grep -i '/lat_lon=' "${tmpfile}" | head -n1 | sed -E 's/.*lat_lon="?([^"]+)"?.*/\1/' | tr -d '\n' || echo "")
    isolate=$(grep -i '/isolate=' "${tmpfile}" | head -n1 | sed -E 's/.*isolate="?([^"]+)"?.*/\1/' | tr -d '\n' || echo "N/A")
    strain=$(grep -i '/strain=' "${tmpfile}" | head -n1 | sed -E 's/.*strain="?([^"]+)"?.*/\1/' | tr -d '\n' || echo "N/A")
    cultivar=$(grep -i '/cultivar=' "${tmpfile}" | head -n1 | sed -E 's/.*cultivar="?([^"]+)"?.*/\1/' | tr -d '\n' || echo "N/A")
    collection_date=$(grep -i '/collection_date=' "${tmpfile}" | head -n1 | sed -E 's/.*collection_date="?([^"]+)"?.*/\1/' | tr -d '\n' || echo "N/A")
    host=$(grep -i '/host=' "${tmpfile}" | head -n1 | sed -E 's/.*host="?([^"]+)"?.*/\1/' | tr -d '\n' || echo "N/A")
    tissue=$(grep -i '/tissue_type=' "${tmpfile}" | head -n1 | sed -E 's/.*tissue_type="?([^"]+)"?.*/\1/' | tr -d '\n' || echo "N/A")
    
    # Parse coordinates if found in GenBank
    if [[ -n "$lat_lon" ]]; then
        read latitude longitude < <(parse_latlon "$lat_lon")
        echo "    ✓ [Method 1: GenBank] Coordinates: $latitude, $longitude" >&2
    fi
    
    # ========================================================================
    # METHOD 2: BioSample via NT Accession Link
    # ========================================================================
    echo "  [NT] Linking to BioSample..." >&2
    biosample=$(timeout 17 bash -c "esearch -db nuccore -query '${acc}' 2>/dev/null | \
                elink -target biosample 2>/dev/null | \
                efetch -format docsum 2>/dev/null | \
                xtract -pattern DocumentSummary -element Accession 2>/dev/null | head -n1" || echo "N/A")
    
    if [[ "$biosample" != "N/A" && -n "$biosample" ]]; then
        echo "    ✓ BioSample found: $biosample" >&2
        
        local bs_data=$(timeout 17 efetch -db biosample -id "$biosample" -format docsum 2>/dev/null || echo "")
        
        if [[ -n "$bs_data" ]]; then
            # Try to get lat_lon from BioSample if not already found
            if [[ "$latitude" == "N/A" || "$longitude" == "N/A" ]]; then
                local bs_latlon=$(echo "$bs_data" | xtract -pattern DocumentSummary -block Attribute \
                    -if Attribute@harmonized_name -equals "lat_lon" -element Attribute 2>/dev/null || echo "")
                
                if [[ -n "$bs_latlon" ]]; then
                    read latitude longitude < <(parse_latlon "$bs_latlon")
                    echo "    ✓ [Method 2: BioSample] Coordinates: $latitude, $longitude" >&2
                fi
            fi
            
            # Fill in missing fields from BioSample
            if [[ "$geo_loc" == "N/A" ]]; then
                geo_loc=$(echo "$bs_data" | xtract -pattern DocumentSummary -block Attribute \
                    -if Attribute@harmonized_name -equals "geo_loc_name" -element Attribute 2>/dev/null || echo "N/A")
            fi
            
            if [[ "$country" == "N/A" ]]; then
                country=$(echo "$bs_data" | xtract -pattern DocumentSummary -block Attribute \
                    -if Attribute@attribute_name -equals "country" -element Attribute 2>/dev/null || echo "N/A")
            fi
            
            if [[ "$isolate" == "N/A" ]]; then
                isolate=$(echo "$bs_data" | xtract -pattern DocumentSummary -block Attribute \
                    -if Attribute@attribute_name -equals "isolate" -element Attribute 2>/dev/null || echo "N/A")
            fi
            
            if [[ "$strain" == "N/A" ]]; then
                strain=$(echo "$bs_data" | xtract -pattern DocumentSummary -block Attribute \
                    -if Attribute@attribute_name -equals "strain" -element Attribute 2>/dev/null || echo "N/A")
            fi
            
            if [[ "$collection_date" == "N/A" ]]; then
                collection_date=$(echo "$bs_data" | xtract -pattern DocumentSummary -block Attribute \
                    -if Attribute@attribute_name -equals "collection_date" -element Attribute 2>/dev/null || echo "N/A")
            fi
            
            if [[ "$host" == "N/A" ]]; then
                host=$(echo "$bs_data" | xtract -pattern DocumentSummary -block Attribute \
                    -if Attribute@attribute_name -equals "host" -element Attribute 2>/dev/null || echo "N/A")
            fi
        fi
    fi
    
    # Clean up fields
    [[ -z "$organism" ]] && organism="N/A"
    [[ -z "$biosample" ]] && biosample="N/A"
    [[ -z "$country" ]] && country="N/A"
    [[ -z "$geo_loc" ]] && geo_loc="N/A"
    [[ -z "$isolate" ]] && isolate="N/A"
    [[ -z "$strain" ]] && strain="N/A"
    [[ -z "$cultivar" ]] && cultivar="N/A"
    [[ -z "$collection_date" ]] && collection_date="N/A"
    [[ -z "$host" ]] && host="N/A"
    [[ -z "$tissue" ]] && tissue="N/A"
    
    # Output: BioSample, Organism, Country, Geo_Location, Latitude, Longitude, Isolate, Strain, Cultivar, Collection_Date, Host, Tissue, Platform, Library
    echo -e "${biosample}\t${organism}\t${country}\t${geo_loc}\t${latitude}\t${longitude}\t${isolate}\t${strain}\t${cultivar}\t${collection_date}\t${host}\t${tissue}\t${platform}\t${library}"
}

################################################################################
# SRA Accession Extraction (Comprehensive Multi-Method)
################################################################################

extract_sra_metadata() {
    local acc="$1"
    
    echo "  [SRA] Fetching run information..." >&2
    
    # Initialize all variables
    local organism="N/A"
    local biosample="N/A"
    local country="N/A"
    local geo_loc="N/A"
    local latitude="N/A"
    local longitude="N/A"
    local isolate="N/A"
    local strain="N/A"
    local cultivar="N/A"
    local collection_date="N/A"
    local host="N/A"
    local tissue="N/A"
    local platform="N/A"
    local library="N/A"
    
    # ========================================================================
    # METHOD 1: SRA RunInfo (Fast method for basic metadata)
    # ========================================================================
    local runinfo=$(timeout 17 bash -c "esearch -db sra -query '$acc' 2>/dev/null | \
                    efetch -format runinfo 2>/dev/null | tail -n1" || echo "")
    
    if [[ -z "$runinfo" ]]; then
        echo "  ✗ Failed to fetch SRA runinfo" >&2
        echo -e "${biosample}\t${organism}\t${country}\t${geo_loc}\t${latitude}\t${longitude}\t${isolate}\t${strain}\t${cultivar}\t${collection_date}\t${host}\t${tissue}\t${platform}\t${library}"
        return 1
    fi
    
    # Parse CSV runinfo
    biosample=$(echo "$runinfo" | cut -d',' -f26)
    organism=$(echo "$runinfo" | cut -d',' -f29)
    platform=$(echo "$runinfo" | cut -d',' -f19)
    library=$(echo "$runinfo" | cut -d',' -f13)
    
    echo "    ✓ [Method 1: RunInfo] Platform: $platform, Library: $library" >&2
    
    # ========================================================================
    # METHOD 2: BioSample via SRA Accession (for geographic data)
    # ========================================================================
    if [[ -n "$biosample" && "$biosample" != "N/A" ]]; then
        echo "  [SRA] Fetching BioSample: $biosample" >&2
        
        local bs_data=$(timeout 17 efetch -db biosample -id "$biosample" -format docsum 2>/dev/null || echo "")
        
        if [[ -n "$bs_data" ]]; then
            # Extract geographic location using harmonized_name
            geo_loc=$(echo "$bs_data" | xtract -pattern DocumentSummary -block Attribute \
                -if Attribute@harmonized_name -equals "geo_loc_name" -element Attribute 2>/dev/null || echo "N/A")
            
            # Extract lat_lon using harmonized_name (PREFERRED METHOD for SRA)
            local bs_latlon=$(echo "$bs_data" | xtract -pattern DocumentSummary -block Attribute \
                -if Attribute@harmonized_name -equals "lat_lon" -element Attribute 2>/dev/null || echo "")
            
            if [[ -n "$bs_latlon" ]]; then
                read latitude longitude < <(parse_latlon "$bs_latlon")
                echo "    ✓ [Method 2: BioSample harmonized] Coordinates: $latitude, $longitude" >&2
            fi
            
            # Extract other attributes
            country=$(echo "$bs_data" | xtract -pattern DocumentSummary -block Attribute \
                -if Attribute@attribute_name -equals "country" -element Attribute 2>/dev/null || echo "N/A")
            
            isolate=$(echo "$bs_data" | xtract -pattern DocumentSummary -block Attribute \
                -if Attribute@attribute_name -equals "isolate" -element Attribute 2>/dev/null || echo "N/A")
            
            strain=$(echo "$bs_data" | xtract -pattern DocumentSummary -block Attribute \
                -if Attribute@attribute_name -equals "strain" -element Attribute 2>/dev/null || echo "N/A")
            
            cultivar=$(echo "$bs_data" | xtract -pattern DocumentSummary -block Attribute \
                -if Attribute@attribute_name -equals "cultivar" -element Attribute 2>/dev/null || echo "N/A")
            
            collection_date=$(echo "$bs_data" | xtract -pattern DocumentSummary -block Attribute \
                -if Attribute@attribute_name -equals "collection_date" -element Attribute 2>/dev/null || echo "N/A")
            
            host=$(echo "$bs_data" | xtract -pattern DocumentSummary -block Attribute \
                -if Attribute@attribute_name -equals "host" -element Attribute 2>/dev/null || echo "N/A")
            
            tissue=$(echo "$bs_data" | xtract -pattern DocumentSummary -block Attribute \
                -if Attribute@attribute_name -equals "tissue" -element Attribute 2>/dev/null || echo "N/A")
        fi
    fi
    
    # Clean up fields
    [[ -z "$biosample" ]] && biosample="N/A"
    [[ -z "$organism" ]] && organism="N/A"
    [[ -z "$platform" ]] && platform="N/A"
    [[ -z "$library" ]] && library="N/A"
    [[ -z "$country" ]] && country="N/A"
    [[ -z "$geo_loc" ]] && geo_loc="N/A"
    [[ -z "$isolate" ]] && isolate="N/A"
    [[ -z "$strain" ]] && strain="N/A"
    [[ -z "$cultivar" ]] && cultivar="N/A"
    [[ -z "$collection_date" ]] && collection_date="N/A"
    [[ -z "$host" ]] && host="N/A"
    [[ -z "$tissue" ]] && tissue="N/A"
    
    # Output: BioSample, Organism, Country, Geo_Location, Latitude, Longitude, Isolate, Strain, Cultivar, Collection_Date, Host, Tissue, Platform, Library
    echo -e "${biosample}\t${organism}\t${country}\t${geo_loc}\t${latitude}\t${longitude}\t${isolate}\t${strain}\t${cultivar}\t${collection_date}\t${host}\t${tissue}\t${platform}\t${library}"
}

################################################################################
# Main Execution
################################################################################

echo "================================================================================" >&2
echo "UNIFIED LOCATION AND METADATA EXTRACTION" >&2
echo "================================================================================" >&2
echo "Start Time: $(date)" >&2
echo "Working Directory: $(pwd)" >&2
echo "" >&2

# Validate required tools
echo "Validating required tools..." >&2
for tool in esearch elink efetch xtract; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "Error: Required tool '$tool' not found. Please install NCBI EDirect utilities." >&2
        exit 1
    fi
done
echo "✓ All required tools available" >&2
echo "" >&2

# Create output file with comprehensive header
cat > "${FINAL_OUTPUT}" << 'EOF'
Accession	Source	BioSample	Organism	Country	Geo_Location	Latitude	Longitude	Isolate	Strain	Cultivar	Collection_Date	Host	Tissue	Platform	Library
EOF

total_accs=0
processed=0
with_coords=0
nt_count=0
sra_count=0

# Count total accessions
if [[ -f "${NT_ACCESSIONS}" && -s "${NT_ACCESSIONS}" ]]; then
    nt_count=$(wc -l < "${NT_ACCESSIONS}")
    total_accs=$((total_accs + nt_count))
fi

if [[ -f "${SRA_ACCESSIONS}" && -s "${SRA_ACCESSIONS}" ]]; then
    sra_count=$(wc -l < "${SRA_ACCESSIONS}")
    total_accs=$((total_accs + sra_count))
fi

if [[ $total_accs -eq 0 ]]; then
    echo "Error: No accessions found in input files." >&2
    echo "  NT accessions file: ${NT_ACCESSIONS}" >&2
    echo "  SRA accessions file: ${SRA_ACCESSIONS}" >&2
    exit 1
fi

echo "Total accessions to process: ${total_accs}" >&2
echo "  NT accessions:  ${nt_count}" >&2
echo "  SRA accessions: ${sra_count}" >&2
echo "" >&2

################################################################################
# Process NT Accessions
################################################################################

if [[ -f "${NT_ACCESSIONS}" && -s "${NT_ACCESSIONS}" ]]; then
    echo "================================================================================" >&2
    echo "PROCESSING NT ACCESSIONS (${nt_count} total)" >&2
    echo "================================================================================" >&2
    echo "" >&2
    
    while IFS= read -r acc || [[ -n "$acc" ]]; do
        # Skip empty lines and comments
        [[ -z "$acc" || "$acc" =~ ^# ]] && continue
        
        # Clean accession (remove whitespace)
        acc=$(echo "$acc" | tr -d '[:space:]')
        
        # Skip if still empty after cleaning
        [[ -z "$acc" ]] && continue
        
        processed=$((processed + 1))
        echo "--------------------------------------------------------------------------------" >&2
        echo "[${processed}/${total_accs}] Processing NT Accession: ${acc}" >&2
        echo "--------------------------------------------------------------------------------" >&2
        
        # Extract metadata (don't let failures stop the loop)
        if metadata=$(extract_nt_metadata "$acc" 2>&1); then
            # Check if has coordinates
            if echo "$metadata" | cut -f5 | grep -qv "N/A"; then
                with_coords=$((with_coords + 1))
            fi
            # Write to output
            echo -e "${acc}\tNT\t${metadata}" >> "${FINAL_OUTPUT}"
        else
            # Write N/A row if extraction failed
            echo -e "${acc}\tNT\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A" >> "${FINAL_OUTPUT}"
            echo "  ✗ Failed to extract metadata, wrote N/A row" >&2
        fi
        
        # Rate limiting to avoid NCBI throttling
        sleep 0.5
        echo "" >&2
        
    done < "${NT_ACCESSIONS}"
fi

################################################################################
# Process SRA Accessions
################################################################################

if [[ -f "${SRA_ACCESSIONS}" && -s "${SRA_ACCESSIONS}" ]]; then
    echo "================================================================================" >&2
    echo "PROCESSING SRA ACCESSIONS (${sra_count} total)" >&2
    echo "================================================================================" >&2
    echo "" >&2
    
    while IFS= read -r acc || [[ -n "$acc" ]]; do
        # Skip empty lines and comments
        [[ -z "$acc" || "$acc" =~ ^# ]] && continue
        
        # Clean accession (remove whitespace)
        acc=$(echo "$acc" | tr -d '[:space:]')
        
        # Skip if still empty after cleaning
        [[ -z "$acc" ]] && continue
        
        processed=$((processed + 1))
        echo "--------------------------------------------------------------------------------" >&2
        echo "[${processed}/${total_accs}] Processing SRA Accession: ${acc}" >&2
        echo "--------------------------------------------------------------------------------" >&2
        
        # Extract metadata (don't let failures stop the loop)
        if metadata=$(extract_sra_metadata "$acc" 2>&1); then
            # Check if has coordinates
            if echo "$metadata" | cut -f5 | grep -qv "N/A"; then
                with_coords=$((with_coords + 1))
            fi
            # Write to output
            echo -e "${acc}\tSRA\t${metadata}" >> "${FINAL_OUTPUT}"
        else
            # Write N/A row if extraction failed
            echo -e "${acc}\tSRA\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A" >> "${FINAL_OUTPUT}"
            echo "  ✗ Failed to extract metadata, wrote N/A row" >&2
        fi
        
        # Rate limiting to avoid NCBI throttling
        sleep 0.5
        echo "" >&2
        
    done < "${SRA_ACCESSIONS}"
fi

################################################################################
# Summary and Statistics
################################################################################

echo "================================================================================" >&2
echo "EXTRACTION COMPLETE!" >&2
echo "================================================================================" >&2
echo "End Time: $(date)" >&2
echo "" >&2
echo "Summary Statistics:" >&2
echo "  Total accessions processed:        ${processed}" >&2
echo "  NT accessions processed:           ${nt_count}" >&2
echo "  SRA accessions processed:          ${sra_count}" >&2
echo "  Accessions with coordinates:       ${with_coords}" >&2
echo "  Accessions without coordinates:    $((processed - with_coords))" >&2
if [[ $processed -gt 0 ]]; then
    echo "  Coordinate success rate:           $(awk "BEGIN {printf \"%.1f%%\", ($with_coords/$processed)*100}")" >&2
fi
echo "" >&2
echo "Output File: ${FINAL_OUTPUT}" >&2
echo "Log File: ${LOG_FILE}" >&2
echo "" >&2

# Display sample results
echo "================================================================================" >&2
echo "SAMPLE RESULTS WITH COORDINATES" >&2
echo "================================================================================" >&2
head -1 "${FINAL_OUTPUT}" >&2
echo "--------------------------------------------------------------------------------" >&2
awk -F'\t' '$7 != "N/A" && NR > 1' "${FINAL_OUTPUT}" | head -10 >&2

echo "" >&2
echo "================================================================================" >&2
echo "SAMPLE RESULTS WITHOUT COORDINATES (showing other metadata)" >&2
echo "================================================================================" >&2
head -1 "${FINAL_OUTPUT}" >&2
echo "--------------------------------------------------------------------------------" >&2
awk -F'\t' '$7 == "N/A" && NR > 1' "${FINAL_OUTPUT}" | head -5 >&2

echo "" >&2
echo "================================================================================" >&2
echo "Done! Check the output file for complete results." >&2
echo "================================================================================" >&2
