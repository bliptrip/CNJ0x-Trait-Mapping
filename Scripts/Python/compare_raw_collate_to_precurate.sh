#!/bin/sh

DATA_FOLDER_PREFIX=$1
DERIVED_FOLDER=$2
raw_collate=$3
curate=$4
output_template=$5
output_file=$6
curate_row_start=$7

#First copy the template the specified output file.
raw_collate="${DATA_FOLDER_PREFIX}/${DERIVED_FOLDER}/${raw_collate}"
curate="${DATA_FOLDER_PREFIX}/${DERIVED_FOLDER}/${curate}"
output_template="${DATA_FOLDER_PREFIX}/${DERIVED_FOLDER}/${output_template}"
output_file="${DATA_FOLDER_PREFIX}/${DERIVED_FOLDER}/${output_file}"

cp -f "${output_template}" "${output_file}"
if [ $# -eq 7 ]; then
    cmd="./compare_raw_collate_to_precurate.py -r \"${raw_collate}\" --raw_ws \"Combined CNJ Upright Data\" -c \"${curate}\" --curated_ws \"Combined CNJ Upright Data\" -o \"${output_file}\" --output_ws \"Combined CNJ Upright Data\" --curate_row_start=$curate_row_start --convert_year"
else
    cmd="./compare_raw_collate_to_precurate.py -r \"${raw_collate}\" --raw_ws \"Combined CNJ Upright Data\" -c \"${curate}\" --curated_ws \"Combined CNJ Upright Data\" -o \"${output_file}\" --output_ws \"Combined CNJ Upright Data\" --curate_row_start=$curate_row_start"
fi 
echo $cmd
eval $cmd 2> errors
