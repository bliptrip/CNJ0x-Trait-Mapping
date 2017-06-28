#!/bin/sh

DATA_FOLDER_PREFIX=$1
raw_collate=$2
curate=$3
output_template=$4
output_file=$5

#First copy the template the specified output file.
raw_collate="${DATA_FOLDER_PREFIX}/Compiled/${raw_collate}"
curate="${DATA_FOLDER_PREFIX}/Compiled/${curate}"
output_template="${DATA_FOLDER_PREFIX}/Compiled/${output_template}"
output_file="${DATA_FOLDER_PREFIX}/Compiled/${output_file}"

cp -f "${output_template}" "${output_file}"
cmd="./compare_raw_collate_to_precurate.py -r \"${raw_collate}\" --raw_ws \"Combined CNJ Upright Data\" -c \"${curate}\" --curated_ws \"Combined CNJ Upright Data\" -o \"${output_file}\" --output_ws \"Combined CNJ Upright Data\""
echo $cmd
eval $cmd 2> errors
