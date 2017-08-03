#!/bin/sh

#Fixes all accessions in a single bash script
echo "This script is supposed to be run from the Script folder in order to function correctly with paths."

DATA_FOLDER_PREFIX=$1
RAW_DATA_PREFIX=$2
DIFF_FILE_SUFFIX=$3

IFS=','; while read file commit; do
    #Step 1: Pull the repository original version of file from git
    FULL_FILENAME=${DATA_FOLDER_PREFIX}/${RAW_DATA_PREFIX}/${file}
    FULL_FILENAME_BASE=${FULL_FILENAME%.*}
    ORIG_FILENAME=${FULL_FILENAME_BASE}.${DIFF_FILE_SUFFIX}.xlsx
    cmd="git checkout $commit -- \"$FULL_FILENAME\""
    echo $cmd
    eval $cmd
    cmd="cp \"$FULL_FILENAME\" \"$ORIG_FILENAME\""
    echo $cmd
    eval $cmd
    cmd="git checkout HEAD -- \"$FULL_FILENAME\""
    echo $cmd
    eval $cmd
    cmd="./diff_spreadsheets.py -1 \"$ORIG_FILENAME\" --file_1_ws \"Sheet1\" -2 \"$FULL_FILENAME\" --file_2_ws \"Sheet1\" -o \"${FULL_FILENAME_BASE}.diffs.xlsx\" --output_ws \"Sheet1\""
    echo $cmd
    eval $cmd
    cmd="rm -f \"$ORIG_FILENAME\""
    echo $cmd
    eval $cmd
done < diff_all_spreadsheets.txt
