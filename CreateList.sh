#!/bin/bash
# Script creates an array for the later programs.
NUMBOFFASTA=$(ls *.faa | wc -l)
NUMBOFGENOME=$(ls *.fa | wc -l)
TOTAL=$(($NUMBOFFASTA*$NUMBOFGENOME))
thiscounterofstuff=1
# This section is for generating a it a job array
# echo source /etc/profile.d/modules.sh>>Masterlist.sh
# echo module load apps/repeatmasker/4.0.5/gcc-4.4.7+trf-4.07b+rmblast-2.2.27>>Masterlist.sh
# echo module load apps/ncbiblast/2.2.29/gcc-4.4.7>>Masterlist.sh
# rm Masterlist.sh
# for fastafile in *.faa; do
# 	for genomedb in *fa; do
# 		echo sh Probe_Check.sh ${genomedb%%.*} $fastafile>>Masterlist.sh
# 	done
# done
rm Masterlist.txt
for fastafile in *.faa; do
	for genomedb in *fa; do
		echo ${genomedb%%.*} $fastafile>>Masterlist.txt
	done
done
