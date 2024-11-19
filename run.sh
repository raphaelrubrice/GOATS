#!/bin/bash
script_dir=$(dirname "$(realpath "${BASH_SOURCE[0]}")")

show_help() {
    echo "Usage: $(basename "$0") [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -g, --geneset       Path to geneset file. ONE GENE SYMBOL PER LINE. (required)"
    echo "  -obo, --obo         Path to the .obo ontology file (required)"
    echo "  -gaf, --gaf         Path to the .gaf annotation file (required)"
    echo "  -c, --collection    Collection(s) to focus on -> must be part of the following 'MF'/'BP'/'CC'. e.g 'MF,CC' or 'MF' (default: 'MF,BP,CC')"       
    echo "  -threshold, --threshold Threshold value for filtering (default: 0.5)"
    echo "  -sd, --show_descriptions    Show descriptions for GO terms (default: 0)"
    echo "  -save, --save          Save the outputs to a file (default: 1)"
    echo "  -name,--name    Name of the folder created if save equals 1"
    echo "  -h, --help             Show this help message and exit"
    echo ""
    echo "Examples:"
    echo "  $(basename "$0") -g geneset.txt -obo go-basic.obo -gaf goa_human.gaf -threshold 0.7"
    echo "  $(basename "$0") --geneset geneset.txt --obo go-basic.obo --gaf goa_human.gaf --show_descriptions 1"
    echo ""
    echo "Note: Ensure that all required files are provided. Missing arguments will result in an error."
    exit 1
}

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --geneset|-g) geneset="$2"; shift ;;
        --obo|-obo) obo_file="$2"; shift ;;
        --gaf|-gaf) gaf_file="$2"; shift ;;
        -c|--collection) collection="$2"; shift ;;
        --threshold|-threshold) threshold="$2"; shift ;;
        --show_descriptions|-sd) show_desc="$2"; shift ;;
        --save|-save) save="$2"; shift ;;
        -name|--name) name="$2"; shift ;;
        --help|-h) show_help ;;
        *) echo "Unknown parameter passed" $1; show_help ;;
    esac
    shift
done

geneset="${geneset:-0}"
obo_file="${obo_file:-0}"
gaf_file="${gaf_file:-0}"
collection="${collection:-"MF,BP,CC"}"
threshold="${threshold:-0.5}"
show_desc="${show_desc:-0}"
save="${save:-1}"
name="${name:-None}"

if [[ -z $obo_file ]] || [[ -z $gaf_file ]] || [[ -z $geneset ]]; then
    echo "Error, you either did not specify the path to the .OBO or .GAF files; or you did not specify any gene"
    show_help
fi

current_date=$(date '+%Y-%m-%d') # Format: YYYY-MM-DD
current_time=$(date '+%H:%M:%S') # Format: HH:MM:SS
log_file="$script_dir/log_${current_date}_${current_time}.txt"
cat > "$log_file" <<EOL
GENESET=$geneset
OBO_FILE=$obo_file
GAF_FILE=$gaf_file
COLLECTIONS=$collection
THRESHOLD=$threshold
SHOW_DESC=$show_desc
SAVE=$save
FOLDERNAME=$name
TIMESTAMP="${current_date}_${current_time}"
EOL

config_file="$script_dir/config.txt"
cat > "$config_file" <<EOL
GENESET=$geneset
OBO_FILE=$obo_file
GAF_FILE=$gaf_file
COLLECTIONS=$collection
THRESHOLD=$threshold
SHOW_DESC=$show_desc
SAVE=$save
FOLDERNAME=$name
TIMESTAMP="${current_date}_${current_time}"
EOL

python=$(command -v python)
$python "$script_dir/goats.py" "$config_file"