#!/bin/bash

extension=${1##*.} 
basename=${1%.*}

echo $#

if [ $# -ge 3 ]; then
    let POINTSIZE=$3
else
    let POINTSIZE=75
fi

if [ $# -ge 4 ]; then
    let OFFSET=$4
else
    let OFFSET=25
fi

text_offset="$OFFSET,$((OFFSET+((POINTSIZE*2)/3)))"

cmd="convert -bordercolor white -border 2 $1 ${basename}.$2.tmp.${extension}"
echo $cmd
eval $cmd
cmd="convert -fill black -font Roboto -weight Bold -pointsize $POINTSIZE -draw \"text $text_offset '$2' gravity 'Center'\" -bordercolor black -border 2 ${basename}.$2.tmp.${extension} ${basename}.$2.${extension}"
echo $cmd
eval $cmd
rm ${basename}.$2.tmp.${extension} 
