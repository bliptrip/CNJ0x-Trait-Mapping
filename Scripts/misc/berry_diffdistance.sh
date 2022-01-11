#!/bin/bash

shapes=$*
cmd="./misc/berry_diffdistance.py -p ../Data/phenotypic\ data/DerivedData/berry_templates -o ../Data/publication/figures/ $shapes"
echo $cmd
eval $cmd

