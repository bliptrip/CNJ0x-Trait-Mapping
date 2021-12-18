#!/bin/sh

./annotatefig1.sh jpegs/IMG_2432.jpg A
./annotatefig1.sh jpegs/IMG_2432_delineated.jpg B
./annotatefig1.sh jpegs/IMG_2427.jpg C
./annotatefig1.sh jpegs/IMG_2420.jpg D
montage -background white -geometry 4032x3024+30+15 -bordercolor black -border 5 -tile 2x2 -label "A" jpegs/IMG_2432.A.jpg jpegs/IMG_2432_delineated.B.jpg jpegs/IMG_2427.C.jpg jpegs/IMG_2420.D.jpg "Fig. 1. Phenotypes.png"
