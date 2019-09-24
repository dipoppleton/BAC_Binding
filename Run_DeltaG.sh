#!/bin/bash
#$ -o test_RM.$JOB_ID.output
#$ -j y
#$ -N Hybrid_Calc
#$ -t 1-1000000
#$ -cwd
#$ -m ae
#$ -e $JOB_NAME.$JOB_ID."$TASK_ID"_error.txt 
#$ -o $JOB_NAME.$JOB_ID."$TASK_ID"_output.txt
#$ -q byslot.q@node08,byslot.q@node03,byslot.q@node04
source /etc/profile.d/modules.sh
module load apps/repeatmasker/4.0.5/gcc-4.4.7+trf-4.07b+rmblast-2.2.27
module load apps/ncbiblast/2.2.29/gcc-4.4.7
i=$SGE_TASK_ID
P=`awk "NR==$i" Masterlist.txt`
sh Probe_Check.sh $P 





