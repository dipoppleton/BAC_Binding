#!/usr/bin/env sh

for fastafile in *.faa; do
	echo "Genome,DeltaG,RPDG_5000,RPDG_1000,RPDG_500,RPDG_0," >${fastafile%%.*}_results.csv
# for when repeat mask is used 	echo "Genome,DeltaG,DeltaG_Repeats,RPDG_5000,RPDG_1000,RPDG_500,RPDG_0,RPDGR_5000,RPDGR_1000,RPDGR_500,RPDGR_0" >${fastafile%%.*}_results.csv
	for genomedb in *.fa; do
		FOLDER=${genomedb%%.*}_${fastafile%%.*}
		cd $FOLDER
		cd RESULTS
		printf "${genomedb%%.*}" >>../../${fastafile%%.*}_results.csv
		printf "," >>../../${fastafile%%.*}_results.csv
		cat ${fastafile%%.*}.out >>../../${fastafile%%.*}_results.csv
		printf "," >>../../${fastafile%%.*}_results.csv
		perl ../../Nonspecific_calc.pl *.NSR.csv >>../../${fastafile%%.*}_results.csv
		echo "">>../../${fastafile%%.*}_results.csv
		cd ..
		cd ..
	done
done

