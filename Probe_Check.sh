#!/usr/bin/env sh

module load apps/repeatmasker/4.0.5/gcc-4.4.7+trf-4.07b+rmblast-2.2.27
module load apps/ncbiblast/2.2.29/gcc-4.4.7


###Start with a genome input and a file input (this program should be looped for each query sequence).
##USAGE Probe_Check.sh Vibrio RNApoltest.faa

### Genome must be present at specified location and blast formated. ncbi BLAST+ must be installed.
GENOME_DB=$1
TESTDB=$GENOME_DB".nsq"
#if [ ! -e $TESTDB ]; then
#	echo "No BLAST database found at "$GENOME_DB
#	echo "Generate it with makeblastdb -dbtype nucl -out FILEINWITHOUTEXTENTION -parse_seqids -in FILEIN.faa "
#	exit
#fi
INPUT_FILE=$2;
if [ ! -e $INPUT_FILE ]; then
	echo "No fasta query file found at "$INPUT_FILE
	exit
fi
#### GET base name for a further output files
BASE_NAME="${INPUT_FILE%.*}"
echo $BASE_NAME
### Make directory to dump data into
#mkdir $BASE_NAME

MAINDIR=$(pwd)
WORKDIR=$MAINDIR"/"$GENOME_DB"_"$BASE_NAME"/"
echo $WORKDIR;
if [ -d $WORKDIR ]; 
then
	rm -R $WORKDIR
fi
mkdir $WORKDIR
cp $INPUT_FILE $WORKDIR
cd $WORKDIR
### query must only contain a single entry or else it will mess up the whole pipeline.
INPUT_FILE=$2;
if [ ! -e $INPUT_FILE ]; then
	echo "No fasta query file found at "$INPUT_FILE
	exit
fi
## Ensure that there is only one sequence in the query
NUMBEROFSEQ=$(grep -c "^>" $INPUT_FILE)
if [ $NUMBEROFSEQ -gt 1 ]; then
	echo "Too many sequences in "$INPUT_FILE
	echo "Exiting"
	exit
fi
if [ $NUMBEROFSEQ -lt 1 ]; then
	echo "No sequences in "$INPUT_FILE
	echo "Exiting"
	exit
fi



#### Perform BLAST
#blastn
blastn -db $MAINDIR"/"$GENOME_DB -query $INPUT_FILE -out $BASE_NAME.bout -evalue 1e-30 -max_target_seqs 250 -outfmt "6 qseqid sseqid qstart qend sstart send evalue nident sstrand"
# OUTPUT queryID Subject ID QueryStart QueryEnd SubjectStart SubjectEnd evalue numberidentical
head -n 250 $BASE_NAME.bout>temp.bout
cat temp.bout>$BASE_NAME.bout
### Determine how many positions were not aligned in the high scoring in the forward direction,
###Get the first position of the genomic in the alignment
STARTOFSUBJ=$(head -n 1 < $BASE_NAME.bout |cut -f5)
###Get the last position of the genomic in the alignment
ENDOFSUBJ=$(head -n 1 < $BASE_NAME.bout |cut -f6)
###Get the subject ID
SEQID=$(head -n 1 < $BASE_NAME.bout |cut -f2)
PROBESTART=$(head -n 1 < $BASE_NAME.bout |cut -f3)
PROBEND=$(head -n 1 < $BASE_NAME.bout |cut -f4)
STRAND=$(head -n 1 < $BASE_NAME.bout |cut -f9)
PROBELENGTH=$(($(wc -m < $INPUT_FILE) - $(head -n 1 < $INPUT_FILE |wc -m) - $(tr -d -c '\n\r' < $INPUT_FILE | wc -c)))
FINALSTART=$STARTOFSUBJ
FINALEND=$ENDOFSUBJ

if [ $STRAND == "minus" ]
then
	TEMPLENGTH=$(($PROBELENGTH-$PROBESTART+PROBEND))
else
	TEMPLENGTH=$(($PROBELENGTH+$PROBESTART-PROBEND))
fi


while IFS="" read -r p || [ -n "$p" ]
do

	##Split the line into variables
	STARTOFOTHER=$(echo $p | awk '{print $5}')
	ENDOFOTHER=$(echo $p | awk '{print $6}')
	SEQIDOTHER=$(echo $p | awk '{print $2}')
	OTHERSTRAND=$(echo $p | awk '{print $9}')



	if [ "$SEQID" = "$SEQIDOTHER" ]
	then
	#	echo " seqid are equal"
		if [ $OTHERSTRAND != $STRAND ]
		then
			TEST="HERE"
# 			echo "EXCLUDE: Different direction same id"
		elif [ $STRAND == "minus" ]
		then
			TEST="HERE"		
			MAXSTART=$((FINALSTART+$TEMPLENGTH))
			MAXEND=$((FINALEND-$TEMPLENGTH))
# 		echo "MAXSTART  "$MAXSTART
# 		echo "MAXEND    "$MAXEND
		#	echo "plus strand WOOT"
		
		####ADD GOES OVER REGION HERE
			if [ $STARTOFOTHER -le $FINALSTART ] && [ $STARTOFOTHER -ge $FINALEND ]
			then
					####This section goes for if it startswithin the sequence
# 				echo "Starts within range so keep the start the same"	
				if [ $FINALEND -lt $ENDOFOTHER ] 
				then
					TEST="HERE"				
# 					echo "EXCLUDED: DUE TO FALLING WITHIN THE SEQUENCe"
				else
# 					echo "INCLUDED: SAME START, same or DIFFERENT END"
					FINALEND=$ENDOFOTHER				
				fi				
			elif [  $ENDOFOTHER -ge $FINALEND  ] && [ $ENDOFOTHER -le $FINALSTART ] 
			then
					####This section goes for if it ends within the sequence
# 				echo "ends within range so keep the end the same"
				if [ $FINALSTART -gt $STARTOFOTHER ] 
				then
					TEST="HERE"
# 					echo "EXCLUDED: DUE TO FALLING WITHIN THE SEQUENCe"
				else
# 					echo "INCLUDED: New or same START, same END"
					FINALSTART=$STARTOFOTHER				
				fi		
			elif [  $MAXSTART -ge $STARTOFOTHER  ] && [ $STARTOFOTHER -gt $FINALSTART ] 
			then
				FINALSTART=$STARTOFOTHER
# 				echo "Starts within reasonable distance"				
			elif [  $MAXEND -le $ENDOFOTHER  ] && [ $ENDOFOTHER -lt $FINALEND ] 
			then
				FINALEND=$ENDOFOTHER				
# 				echo "ends within reasonable distance"
			else
				TEST="HERE"			
			fi
# 		echo "Start  "$FINALSTART
# 		echo "END    "$FINALEND
		TEMPLENGTH=$(($PROBELENGTH-$FINALSTART+$FINALEND))
# 		echo $TEMPLENGTH
		else 
		###Get the swing in either direction.
			MAXSTART=$((FINALSTART-$TEMPLENGTH))
			MAXEND=$((FINALEND+$TEMPLENGTH))
# 		echo "MAXSTART  "$MAXSTART
# 		echo "MAXEND    "$MAXEND
		#	echo "plus strand WOOT"		
			if [ $STARTOFOTHER -ge $FINALSTART ] && [ $STARTOFOTHER -le $FINALEND ]
			then
					####This section goes for if it startswithin the sequence
# 				echo "Starts within range so keep the start the same"	
				if [ $FINALEND -gt $ENDOFOTHER ] 
				then
# 					echo "EXCLUDED: DUE TO FALLING WITHIN THE SEQUENCe"
					TEST="HERE"
				else
# 					echo "INCLUDED: SAME START, same or DIFFERENT END"
					FINALEND=$ENDOFOTHER				
				fi				
			elif [  $ENDOFOTHER -le $FINALEND  ] && [ $ENDOFOTHER -ge $FINALSTART ] 
			then
					####This section goes for if it ends within the sequence
# 				echo "ends within range so keep the end the same"
				if [ $FINALSTART -lt $STARTOFOTHER ] 
				then
# 					echo "EXCLUDED: DUE TO FALLING WITHIN THE SEQUENCe"
					TEST="HERE"
				else
# 					echo "INCLUDED: New or same START, same END"
					FINALSTART=$STARTOFOTHER				
				fi		
			elif [  $MAXSTART -le $STARTOFOTHER  ] && [ $STARTOFOTHER -lt $FINALSTART ] 
			then
				FINALSTART=$STARTOFOTHER
# 				echo "Starts within reasonable distance"				
			elif [  $MAXEND -ge $ENDOFOTHER  ] && [ $ENDOFOTHER -gt $FINALEND ] 
			then
					FINALEND=$ENDOFOTHER				
# 				echo "ends within reasonable distance"
			else
			TEST="HERE"
			fi
# 		echo "Start  "$FINALSTART
# 		echo "END    "$FINALEND
		TEMPLENGTH=$(($PROBELENGTH+$FINALSTART-$FINALEND))
# 		echo $TEMPLENGTH
		fi
	else 
# 	echo "EXCLUDED: Different ID"
			TEST="HERE"
	fi		
done<$BASE_NAME.bout
##############################


if [[ $FINALEND -lt 0 ]]
then
	FINALEND=0
fi

if [[ $FINALSTART -lt 0 ]]
then
	FINALSTART=0
fi

####Extract genomic sequence and 
if [ $FINALSTART -le $FINALEND ]
then
	blastdbcmd -db $MAINDIR"/"$GENOME_DB -entry $SEQID -range $FINALSTART-$FINALEND -out $BASE_NAME".temp.faa" -strand $STRAND
else
	blastdbcmd -db $MAINDIR"/"$GENOME_DB -entry $SEQID -range $FINALEND-$FINALSTART -out $BASE_NAME".temp.faa" -strand $STRAND
fi
cat $INPUT_FILE >> $BASE_NAME".temp.faa"
#### Perform MAFFT alignment
../mafft $BASE_NAME.temp.faa > $BASE_NAME.aligned.faa

##Make the file a nice fasta The best line of code in existence!
awk '$0~/^>/{if(NR>1){print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' "$BASE_NAME.aligned.faa" > "$BASE_NAME.aligned.nice.faa"

## Get the repeated regions
# RepeatMasker -pa 4 -qq $BASE_NAME.temp.faa
# awk '$0~/^>/{if(NR>1){print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' "$BASE_NAME.temp.faa.masked" > "$BASE_NAME.temp.faa.masked.nice.faa"
# perl $MAINDIR"/"repeatconvertor.pl $BASE_NAME.aligned.nice.faa $BASE_NAME.temp.faa.masked.nice.faa $BASE_NAME.rout

rm *faa.out *faa.cat *faa.tbl

rm -r *.temp.*






#### Calculate delta G
perl $MAINDIR"/"Blast2Energy_Calc.pl $BASE_NAME.aligned.nice.faa $BASE_NAME $BASE_NAME.rout DONTPRINT

# rm *.aligned*  
rm *rout

mkdir "RESULTS"
mv *.out *.csv "RESULTS"
cd ..
