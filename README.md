# BAC_Binding
This program calculates a DeltaG for the Binding of a Bacterial Artifical Chromosome (BAC) onto genomic sequences.
This program can take any nucletide sequence and calculate its theoretical binding efficy on any target.
This is usefull for mismatched BACS or using a single probe on divergent species.

Requirments: 
BLAST+ (any version):ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
If including nonspecific binding agents: repeatmasker with RMBlast http://www.repeatmasker.org/ 
mafft :https://mafft.cbrc.jp/alignment/software/ (MAFFT itself must be in the working directory or path)

Outputs:
a folder named Query_target containg the alignments and intermediatory data and a subfolder named results containg a CSV of a position by position score and a text file with the final deltaG value.

Binding energies and variables can be changed within the file Blast2Energy_Calc.pl


Using the program:
Step 1: Setting up inputs
Place the query BAC sequences in the folder with the extension .faa (ex. BC9331.fa, CW293401.fa)
Querys can only be a single fasta sequence per file.
Create a blast database for each nucleotide target genome ensure parseids flag is selected.
$ makeblastdb -dbtype nucl -out FILEINWITHOUTEXTENTION -parse_seqids -in FILEIN.faa
Leave the original multifasta target file with an extension of .fa (ex. hg38.fa, ecoli.fa)

Step A2: Running locally or withoug Job array
Comment out all text which is commented within CreateList.sh. Uncomment all previously commented text.
After files are located in a single folder run $ sh CreateList.sh
All steps will be taken care off. Not recomended for large quantities of comparisons. Sutiable for testing.
For one on one sequences, just run sh Probe_Check_4.sh target query

Step 2: Setting up GSE job array.
After files are located in a single folder run $ CreateList.sh

Step 3: Run Job.
$ sh Run_DeltaG.sh
If you want to run this job with repeat masker (say you are using a repeat blocker in your experiment), Uncomment lines 235-7
This step can take seconds to half an hour, depending on your input size and whether you choose to test repeats.

Analysis:
To combine all of the csv files for each query, run sh organizeresults.sh
This will provide a single csv file containing all the position by position deltaG.
Should you only wish to see the cumulutive deltaG's run sh GenerateSummary.sh

