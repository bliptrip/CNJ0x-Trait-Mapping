#!/bin/sh

#Fixes all accessions in a single bash script
echo "This script is supposed to be run from the Script folder in order to function correctly with paths."

output_template="Data-combined-template.xlsx"
output_file="$1"
DATA_FOLDER_PREFIX="../Data/phenotypic data/NJ Data"

#First copy the template the specified output file.
output_template="${DATA_FOLDER_PREFIX}/Compiled/${output_template}"
output_file="${DATA_FOLDER_PREFIX}/Compiled/${output_file}"

cp -f "${output_template}" "${output_file}"

while read file population year; do
    cmd="./collate_separate_excel_datasets.py -i \"${DATA_FOLDER_PREFIX}/${file}\" --input_ws \"Sheet1\" -o \"${output_file}\" --output_ws \"Combined CNJ Upright Data\" -p $population -y $year"
    echo "$cmd"
    eval "$cmd"
done < collate_datasets.txt
