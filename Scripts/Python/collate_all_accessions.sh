#!/bin/sh

#Fixes all accessions in a single bash script
echo "This script is supposed to be run from the Script folder in order to function correctly with paths."

DATA_FOLDER_PREFIX=$1
RAW_DATA_PREFIX=$2
DERIVED_DATA_PREFIX=$3
output_template=$4
output_file=$5

#First copy the template the specified output file.
output_template="${DATA_FOLDER_PREFIX}/${DERIVED_DATA_PREFIX}/${output_template}"
output_file="${DATA_FOLDER_PREFIX}/${DERIVED_DATA_PREFIX}/${output_file}"

cp -f "${output_template}" "${output_file}"

while read file population year; do
    cmd="./collate_separate_excel_datasets.py -i \"${DATA_FOLDER_PREFIX}/${RAW_DATA_PREFIX}/${file}\" --input_ws \"Sheet1\" -o \"${output_file}\" --output_ws \"Combined CNJ Upright Data\" -p $population -y $year"
    echo "$cmd"
    eval "$cmd"
done < collate_datasets.txt
