#!/usr/bin/env sh
printf "BAC CLONE," > Final_Results.csv
for genomedb in *.fa; do
	printf ${genomedb%%.*}>>Final_Results.csv
	 printf "," >>Final_Results.csv
done
echo "" >>Final_Results.csv
for fastafile in *.faa; do
	printf ${fastafile%%.*}>>Final_Results.csv
	printf "," >>Final_Results.csv
	for genomedb in *.fa; do
		FOLDER=${genomedb%%.*}_${fastafile%%.*}
	#	echo $FOLDER
		if [ ! -d "$FOLDER" ] ; then
			printf 'FAIL,' >>Final_Results.csv
		else
			cd $FOLDER
			if [ ! -d RESULTS ] ; then 
				printf 'FAIL,' >>../Final_Results.csv
	#	       echo "NO"	
			else
				cd RESULTS
	#		echo "Yes"
				awk -F"." '{printf $1}' ${fastafile%%.*}.out>>../../Final_Results.csv
				printf "," >>../../Final_Results.csv
				cd ..
			fi
			cd ..
		fi
	done
		echo "" >>Final_Results.csv
done

