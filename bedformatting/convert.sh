#!/bin/bash

IFS='.' read -ra BED <<< "$1"

sort -V -k1,1 -k2,2n -k3,3n -u $1 | sed 's/\s*$//' | sed 's/\s/,/g' | sed 's/,/\t/'  | sed 's/,/\t/' | sed 's/,/\t/' > ${BED[0]}-sorted.bed

perl bed6_2_bed12.pl -i ${BED[0]}-sorted.bed -o ${BED[0]}-bed12.bed

Rscript --vanilla targets.R ${BED[0]}
