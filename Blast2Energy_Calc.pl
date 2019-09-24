#!/usr/bin/env perl
use Memoize;
memoize('Mismatch_calculate');
$file_in = $ARGV[0];
$file_out = $ARGV[1];
$file_repeat = $ARGV[2];

##this value is whether any print out is displayed

$display = $ARGV[3];
my $Disp_Out = 1;

if ($#ARGV == 3){
	$Disp_Out  = 0;
}
	print "\n"."alignment2Energy_calc.pl $file_in $file_out $file_out \n";

if ($Disp_Out  == 1){
	print "\n".'alignment2Energy_calc.pl INFILE OUtFILE REPEATFILE printoutput?
	this program takes an alignment input and transforms it into an overall delta g for binding and a table of deltags for each position
	File input must be a 2 sequence nucleotide nice FaSta format aka:
	>seq1
	-attattgc
	>seq2
	tattacggc
	*** NIcE FaSta MEaNS NO NEW LINES WItHIN SEQUENcES.
	tto convert do :awk \'$0~/^>/{if(NR>1){print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}\' "filename" 
	This program generates 2 outfiles, 4 if there is no repeat specified all have the basename of OUTFILE with .csv .out or .rp.csv .rp.out
	The CSV outputs a delta G for each position as a comma seperated string of values. (ideal for R), the .out is a single best deltaG for binding
	the third argument is whether any output is printed. Put any string in to annul any std out
	'."\n";
}
###LOaD FILE
open my $handle, '<', $file_in;
chomp(my @sequences = <$handle>);
close $handle;

###
### PLacE SEQUENcES IN StRINgS
foreach my $lines (@sequences) {$lines =~ tr/[A-Z]/[a-z]/;}
@sequence1=split("",$sequences[1]);
chomp @sequence1;
@sequence2=split("",$sequences[3]);
chomp @sequence2;

# if ($Disp_Out  == 1){print "the two sequences are: \n$sequences[0]\n$sequences[2]\n";}


#foreach my $var (@sequence1){print $var." "}
#print "\n";
#foreach my $var (@sequence2){print $var." "}


###########################tHE actUaL LOOKUP.
#####	M=Match	X=Mismatch	G=gap 	 on bottom N=bridged sequence.

###Start by determining whether each pair is a gap, mismatch or match or a human repeat
@deltagbyposition;
@positiondetails=();

my $counter=0;
while($counter<=$#sequence1){
#	print $sequence1[$counter]."_".$sequence2[$counter]."\t";
	#Match
	if ($sequence1[$counter] eq $sequence2[$counter]){
		$positiondetails[$counter]="M";
	#	print "Match\n";
	}elsif($sequence1[$counter] eq "-" or $sequence2[$counter] eq "-"){
	#	print "gap\n";
		$positiondetails[$counter]="G";
	}elsif($sequence2[$counter] eq "N" or $sequence2[$counter] eq "N"){
	#	print "ambigious\n";
		$positiondetails[$counter]="N";
	}else{
	## Mismatch
	#	print "Mismatch\n";
		$positiondetails[$counter]="X";

	}
	$counter++;
}



#### Start calculating binding energies from one end and work your way in. Every time more than one mismatch is found, calculate the total for the gap.

#print "Before calc\n positiondetails : $#positiondetails Sequence1 $#sequence1 Sequence2 $#sequence2 deltagbyposition $#deltagbyposition\n";
$counter=0;
while($counter<$#positiondetails){
#	print $counter."_";
#	print $positiondetails[$counter]."_".$positiondetails[$counter+1]."_";
	if (exists $deltagbyposition[$counter]){
	}elsif ($positiondetails[$counter] eq "N"){
		$deltagbyposition[$counter] = 0;
		#print "foun an N\n";
	
	#check if this is a double mismatch or gap placement.
	}elsif (($positiondetails[$counter] eq "X" or $positiondetails[$counter] eq "G" or $positiondetails[$counter] eq "N")and ($positiondetails[$counter+1] eq "X" or $positiondetails[$counter+1] eq "G" or $positiondetails[$counter+1] eq "N")){
	#Backcheck to an N
#	print "double mismatch or gap\n";
	#	print "foun an XG and it is $positiondetails[$counter]_$positiondetails[$counter+1]\n";
		my $startposition =$counter;
		my $endposition;
		if ($counter == 0){
			$deltagbyposition[$counter] = 0;
		}else{
			my $tempdg;
			my $backcount=0;
			#Count backwards till a match is found, should be minimal.
			until (($backcount > $counter) or ($positiondetails[$counter-$backcount] eq "M")){
				#print "\n";
				if ($positiondetails[$counter-$backcount-1]eq "M"){
					$startposition=($counter-$backcount-1);
				}elsif ($backcount == $counter){
					$deltagbyposition[$counter] = 0;
				}
			
			$backcount++;
			}
			my $forwardcount=0;
			until ((($forwardcount+$counter) > $#positiondetails) or ($positiondetails[$counter+$forwardcount] eq "M")){
				#print "\n";
				if ($positiondetails[$counter+$forwardcount+1]eq "M"){
					$startposition=($counter+$forwardcount+1);
				}elsif ($forwardcount == $#positiondetails){
					$deltagbyposition[$counter] = 0;
				}
			
			$forwardcount++;
			}
			#print "\n";
			if ($deltagbyposition[$counter] ne 0){
			#print "$backcount\t$forwardcount\n";
				my $inputquery = join( '', @sequence1[($counter-$backcount)..($counter+$forwardcount)]);
				my $inputsubject =join( '', @sequence2[($counter-$backcount)..($counter+$forwardcount)]);
				#print "$inputquery\n";
				#print "$inputsubject\n";

				$tempdg=Gap_Calc($inputquery,$inputsubject);
				#print $tempdg."\n\n";
			#print "$tempdg tempdg\n";

				
			}
			#while bla is less than the total do make delta g
			$backcount =$backcount+1; 
			$forwardcount --;
		#	print "\n\n$counter\n";
		#	print "$backcount backcount\n";
		#	print "$forwardcount forwardcount\n";
			my $total=$backcount+$forwardcount-1;
		#	print "$total total\n";
		#	print "$tempdg tempdg\n";

			$tempdg=$tempdg/$total;
		#	print "$tempdg tempdg\n";ß

			while ($backcount<=0){
			#	unless (exists $deltagbyposition[$counter-$backcount]){
					$deltagbyposition[$counter-$backcount]=$tempdg;
		#		}
			#	print $deltagbyposition[$counter-$backcount]."backcount\n";
				$backcount++;
			}
			while ($forwardcount>=0){
				unless (exists $deltagbyposition[$counter+$forwardcount]){
					$deltagbyposition[$counter+$forwardcount]=$tempdg;
				}
				$forwardcount--;
			#	print ($counter+$forwardcount)."forcount\n";

			}
		}
	#Get the correct number of positions


	#check if this is a single gap position.
	
	}elsif ($positiondetails[$counter] eq "G"){
	#	print "single gap\n";

		my $inputquery = join( '', @sequence1[($counter-1)..($counter+1)]);
		my $inputsubject =join( '', @sequence2[($counter-1)..($counter+1)]);
	
		my $tempdg=Gap_Calc($inputquery,$inputsubject);
		#print $tempdg."\n\n";
	$deltagbyposition[$counter]=$tempdg/2;
	$deltagbyposition[$counter-1]=$tempdg/2;

	}else{
#	#print $sequence1[$counter].$sequence1[$counter+1]."_".$sequence2[$counter].$sequence2[$counter+1]."\n";
	##print $positiondetails[$counter].$positiondetails[$counter+1]."_".$positiondetails[$counter].$positiondetails[$counter+1]."\n";
		#Match
		if ($positiondetails[$counter].$positiondetails[$counter+1] eq "MM"){
#print "Match\n";
#			print &Match_calculate($sequence1[$counter], $sequence1[$counter+1])."\n\n";
			$deltagbyposition[$counter]=&Match_calculate($sequence1[$counter], $sequence1[$counter+1]);
			
#		print "Mismatch forward\n";
		}elsif($positiondetails[$counter].$positiondetails[$counter+1] eq "MX"){
		#	print "Mismatch top\n";
			$deltagbyposition[$counter]= &Mismatch_calculate($sequence1[$counter].$sequence1[$counter+1], $sequence2[$counter].$sequence2[$counter+1]);
#		print "Mismatch backward\n";
		}elsif($positiondetails[$counter].$positiondetails[$counter+1] eq "XM"){
#			print "Mismatch Bottom\n";
			$deltagbyposition[$counter]= &Mismatch_calculate($sequence2[$counter+1].$sequence2[$counter], $sequence1[$counter+1].$sequence1[$counter]);
		#	print $deltagbyposition[$counter];
		}elsif($positiondetails[$counter].$positiondetails[$counter+1] eq "XX"){
#		print "Double mismatch\n";
		}else{
#		print "error\n";
		}
	
	

	}
	$counter++;
#	if ($counter >=4070 ){print "$counter _ $#deltagbyposition\n$positiondetails[$counter]\n"};
}
#print "\n";
#print "\n";

# foreach my $var (@deltagbyposition){print $var.","}
# print "\n";
#print "after calc\n positiondetails : $#positiondetails Sequence1 $#sequence1 Sequence2 $#sequence2 deltagbyposition $#deltagbyposition\n";

###Output the list of deltagpositions to outfile
open(FO, '>', $file_out."\.csv");
print FO "Position\,DeltaG\n";
#print "\n";
my $counter = 0;
foreach my $var (@deltagbyposition){
	if ($var eq ""){$var=0}
	#print "$var,";
	print FO "$counter\,$var\n";
	$counter++;
 	}
print FO "\n";
#print "\n";
close FO;
#print "\n";
foreach my $var (@positiondetails){
	#print "$var,";
 }

my $absolutedeltaG=0;
my $postemp1=0;
my $negtemp1=0;
#print "\nValue\tTracks\tabsolute\tpostemp\tnegtem;\n";
my $tracks=0;###Keep track of what we were building on 1=negative  0=positive.(positive is bad)
my $startingone=0;
foreach my $var (@deltagbyposition){

	my $places = 2;
	my $factor = 10**$places;
	$var=int($var * $factor) / $factor;
#Start by not adding anything that is above zero at the start
	if ($startingone==0 and $var >=0){	
	}else{
		$startingone=1;
		if ($tracks==0){
			#check if we are building on a positive
			if ($var >=0){
				#If this is positive just build on it.
				$postemp1=$postemp1+$var;
			}else{
				$tracks=1;
				$negtemp1=$var+$negtemp1;
			}
	
		}elsif($tracks==1){
				#check if we are building on a negative
			if ($var >=0){
				##switch tracks
				$tracks=0;
				if($absolutedeltaG==0){
				#If no value has been set for a benificial interaction, use this one and reset the neg and pos deltaG.
					$absolutedeltaG=$negtemp1;
					$negtemp1=0;
					$postemp1=$var;
				}elsif($negtemp1<$absolutedeltaG){
				##See if the deltaG is larger than the one stored
					if(($postemp1+$absolutedeltaG)>0){
						#If the mismatch between the sequences sucks worse than the original peak, use this one instead.
						$absolutedeltaG=$negtemp1;
						$negtemp1=0;
						$postemp1=$var;
					}else{
						#If not, add all three values togethor and consider it a single peak. Then reset them
						$absolutedeltaG=$negtemp1+$absolutedeltaG+$postemp1;
						$negtemp1=0;
						$postemp1=$var;
					}
				}else{
					##If it is a benificial interaction add it to the absolute
					if (($negtemp1+$postemp1)<0){
						$absolutedeltaG=$negtemp1+$absolutedeltaG+$postemp1;
						$negtemp1=0;
						$postemp1=$var;
					}else{
						$postemp1=$negtemp1+$postemp1+$var;
						$negtemp1=0;
					}
				}
			}else{
				# If it was a negative build on it.
				$negtemp1=$negtemp1+$var;
			}
		}else{
			#print "error";
		}
	}
	#print "$var\t$tracks\t$absolutedeltaG\t$postemp1\t$negtemp1;\n";
}

if($negtemp1<$absolutedeltaG){
	if(($absolutedeltaG+$postemp1)<0){
		$absolutedeltaG=$negtemp1+$absolutedeltaG+$postemp1;
	}else{
		$absolutedeltaG=$negtemp1;
	}
}else{
	if(($negtemp1+$postemp1)<0){
		$absolutedeltaG=$negtemp1+$absolutedeltaG+$postemp1;
	}
}

#print "FINAL IS $absolutedeltaG\n";
open(FO, '>', $file_out."\.out");
print FO "$absolutedeltaG";
#print "\n";
close FO;


###############################THIS SECTION IS BASICLY A REPEAT FOR USING REPEAT MASKER ON


%deltagbyposition=();

if (-e $file_repeat){

	##Get the repeat sequences
	open my $handle, '<', $file_repeat;
	chomp(my @repeats = <$handle>);
	close $handle;
	@repeat1=split(/,/,$repeats[0]);
	@repeat2=split(/,/,$repeats[1]);
	foreach my $var (@repeat1){ 
	#	print "Position=$var\n";
	#	print "Before=$positiondetails[$var]\n";
		unless ($positiondetails[$var] eq "N"){
			$positiondetails[$var]="R";
			$sequence1[$var]="R";
		}
	#	print "After=$positiondetails[$var]\n";
	}
	#print "\n";

	foreach my $var (@repeat2){ 
	#	print "Position=$var\n";
	#	print "Before=$positiondetails[$var]\n";
		unless ($positiondetails[$var] eq "N"){
			$positiondetails[$var]="R";
			$sequence2[$var]="R";

		}
	#	print "After=$positiondetails[$var]\n";
	}
#	print "\n";

foreach my $var (@sequence1){
#print $var."_";
}
#print "\n";
foreach my $var (@sequence2){
#print $var."_";
}
#print "\n";

#print "START LOOKING HERE\n";
######Checked to here!

my $counter=0;
while($counter<$#positiondetails){
	#print $positiondetails[$counter]."_";
	if (exists $repdeltagbyposition[$counter]){
	}elsif ($positiondetails[$counter] eq "N"){
		$repdeltagbyposition[$counter] = 0;
		#print "foun an N\n";
	#}elsif ($positiondetails[$counter] eq "R"){
	#	#print "foun an R\n";
	
	#check if this is a double mismatch or gap placement.
	}elsif (($positiondetails[$counter] eq "X" or $positiondetails[$counter] eq "G" or $positiondetails[$counter] eq "N" or $positiondetails[$counter] eq "R")and ($positiondetails[$counter+1] eq "X" or $positiondetails[$counter+1] eq "G" or $positiondetails[$counter+1] eq "R")){
	#Backcheck to an N
	#	#print "foun an XG and it is $positiondetails[$counter]_$positiondetails[$counter+1]\n";
		my $startposition =$counter;
		my $endposition;
		if ($counter == 0){
			$repdeltagbyposition[$counter] = 0;
		}else{
			my $tempdg;
			my $backcount=0;
			#Count backwards till a match is found, should be minimal.
			until (($backcount > $counter) or ($positiondetails[$counter-$backcount] eq "M")){
				##print "\n";
				if ($positiondetails[$counter-$backcount-1]eq "M"){
					$startposition=($counter-$backcount-1);
				}elsif ($backcount == $counter){
					$repdeltagbyposition[$counter] = 0;
				}
			
			$backcount++;
			}
			my $forwardcount=0;
			until ((($forwardcount+$counter) > $#positiondetails) or ($positiondetails[$counter+$forwardcount] eq "M")){
				##print "\n";
				if ($positiondetails[$counter+$forwardcount+1]eq "M"){
					$startposition=($counter+$forwardcount+1);
				}elsif ($forwardcount == $#positiondetails){
					$repdeltagbyposition[$counter] = 0;
				}
			
			$forwardcount++;
			}
			#print "\n";
			if ($repdeltagbyposition[$counter] ne 0){
			#print "$backcount\t$forwardcount\n";
				my $inputquery = join( '', @sequence1[($counter-$backcount)..($counter+$forwardcount)]);
				my $inputsubject =join( '', @sequence2[($counter-$backcount)..($counter+$forwardcount)]);
				#print "$inputquery\n";
				#print "$inputsubject\n";

				$tempdg=Gap_Calc($inputquery,$inputsubject);
				#print $tempdg."\n\n";
			#print "$tempdg tempdg\n";

				
			}
			#while bla is less than the total do make delta g
			$backcount ++; 
			$forwardcount --;
			my $total=$backcount+$forwardcount-1;
			#print "$total total\n";
			#print "$tempdg tempdg\n";

			$tempdg=$tempdg/$total;
			#print "$tempdg tempdg\n";

			while ($backcount<=0){
				$repdeltagbyposition[$counter-$backcount]=$tempdg;
				$backcount++;
			}
			while ($forwardcount>=0){
				$repdeltagbyposition[$counter+$forwardcount]=$tempdg;
				$forwardcount--;
			}
		}
	#Get the correct number of positions


	#check if this is a single gap position.
	
	}elsif ($positiondetails[$counter] eq "G"){
	
		my $inputquery = join( '', @sequence1[($counter-1)..($counter+1)]);
		my $inputsubject =join( '', @sequence2[($counter-1)..($counter+1)]);
	
		my $tempdg=Gap_Calc($inputquery,$inputsubject);
		#print $tempdg."\n\n";
	$repdeltagbyposition[$counter]=$tempdg/2;
	$repdeltagbyposition[$counter-1]=$tempdg/2;

	}else{
#	#print $sequence1[$counter].$sequence1[$counter+1]."_".$sequence2[$counter].$sequence2[$counter+1]."\n";
	##print $positiondetails[$counter].$positiondetails[$counter+1]."_".$positiondetails[$counter].$positiondetails[$counter+1]."\n";
		#Match
		if ($positiondetails[$counter].$positiondetails[$counter+1] eq "MM"){
	#		print "Match\n";
#			print &Match_calculate($sequence1[$counter], $sequence1[$counter+1])."\n\n";
			$repdeltagbyposition[$counter]=&Match_calculate($sequence1[$counter], $sequence1[$counter+1]);
			
		#Mismatch forward
		}elsif($positiondetails[$counter].$positiondetails[$counter+1] eq "MX"){
		#	print "Mismatch top\n";
			$repdeltagbyposition[$counter]= &Mismatch_calculate($sequence1[$counter].$sequence1[$counter+1], $sequence2[$counter].$sequence2[$counter+1]);
		#Mismatch backward
		}elsif($positiondetails[$counter].$positiondetails[$counter+1] eq "XM"){
	#		print "Mismatch Bottom\n";
			$repdeltagbyposition[$counter]= &Mismatch_calculate($sequence2[$counter+1].$sequence2[$counter], $sequence1[$counter+1].$sequence1[$counter]);
		#	print $repdeltagbyposition[$counter];
		}elsif($positiondetails[$counter].$positiondetails[$counter+1] eq "RM"){
			$repdeltagbyposition[$counter]=0.5;
		}elsif($positiondetails[$counter].$positiondetails[$counter+1] eq "MR"){
			$repdeltagbyposition[$counter]=0.5;
		}else{
#		print "error\n";
		}
	
	

	}
	$counter++;
}






###Output the list of deltagpositions to outfile
open(FO, '>', $file_out."\.rp\.csv");

print FO "Position\,DeltaG\n";
#print "\n";
my $counter = 0;
foreach my $var (@repdeltagbyposition){
	if ($var eq ""){$var=0}
	#print "$var,";
	print FO "$counter\,$var\n";
	$counter++;
 }
print FO "\n";


#print "\n";
close FO;
#print "\n";


$absolutedeltaG=0;
$postemp1=0;
$negtemp1=0;
#print "\nValue\tTracks\tabsolute\tpostemp\tnegtem;\n";
$tracks=0;###Keep track of what we were building on 1=negative  0=positive.(positive is bad)
$startingone=0;
foreach my $var (@repdeltagbyposition){

	my $places = 2;
	my $factor = 10**$places;
	$var=int($var * $factor) / $factor;
#Start by not adding anything that is above zero at the start
	if ($startingone==0 and $var >=0){	
	}else{
		$startingone=1;
		if ($tracks==0){
			#check if we are building on a positive
			if ($var >=0){
				#If this is positive just build on it.
				$postemp1=$postemp1+$var;
			}else{
				$tracks=1;
				$negtemp1=$var+$negtemp1;
			}
	
		}elsif($tracks==1){
				#check if we are building on a negative
			if ($var >=0){
				##switch tracks
				$tracks=0;
				if($absolutedeltaG==0){
				#If no value has been set for a benificial interaction, use this one and reset the neg and pos deltaG.
					$absolutedeltaG=$negtemp1;
					$negtemp1=0;
					$postemp1=$var;
				}elsif($negtemp1<$absolutedeltaG){
				##See if the deltaG is larger than the one stored
					if(($postemp1+$absolutedeltaG)>0){
						#If the mismatch between the sequences sucks worse than the original peak, use this one instead.
						$absolutedeltaG=$negtemp1;
						$negtemp1=0;
						$postemp1=$var;
					}else{
						#If not, add all three values togethor and consider it a single peak. Then reset them
						$absolutedeltaG=$negtemp1+$absolutedeltaG+$postemp1;
						$negtemp1=0;
						$postemp1=$var;
					}
				}else{
					##If it is a benificial interaction add it to the absolute
					if (($negtemp1+$postemp1)<0){
						$absolutedeltaG=$negtemp1+$absolutedeltaG+$postemp1;
						$negtemp1=0;
						$postemp1=$var;
					}else{
						$postemp1=$negtemp1+$postemp1+$var;
						$negtemp1=0;
					}
				}
			}else{
				# If it was a negative build on it.
				$negtemp1=$negtemp1+$var;
			}
		}else{
			#print "error";
		}
	}
	#print "$var\t$tracks\t$absolutedeltaG\t$postemp1\t$negtemp1;\n";
}

if($negtemp1<$absolutedeltaG){
	if(($absolutedeltaG+$postemp1)<0){
		$absolutedeltaG=$negtemp1+$absolutedeltaG+$postemp1;
	}else{
		$absolutedeltaG=$negtemp1;
	}
}else{
	if(($negtemp1+$postemp1)<0){
		$absolutedeltaG=$negtemp1+$absolutedeltaG+$postemp1;
	}
}

#print "FINAL IS $absolutedeltaG\n";
open(FO, '>', $file_out."\.rp\.out");
print FO "$absolutedeltaG";
#print "\n";
close FO;

















}

###############################END REPEAT SECTION




















#memoize(function); ###Speed up a slow fuction if necessary


 # sub Reverese_Compliment{
#  	my $origin_seq =shift;
#  	my $revcomp = reverse $origin_seq;
#  	$revcomp =~ tr/ATGCatgc/TACGtacg/;
#  	return $revcomp;
#  }



sub Mismatch_calculate{
	my $Query = shift;
#	print $Query."query\n";
	my $Subject =shift;
	#print $Subject."Subject\n";
	my $Aligned=substr($Query,0,1);
	my $mismatch=substr($Query,1,1).substr($Subject,1,1);
#	print "KEY is $Aligned \n while mismatch is $mismatch\n";
	if ($Aligned eq "g"){
		my %MX_GAP_G = (
			aa => 'Print "error"',
			at => '0.17',
			ac => '-0.25',
			ag => '0.81',
			ta => '0.45',
			tt => 'Print "error"',
			tc => '-0.59',
			tg => '0.98',
			ca => '0.62',
			ct => '0.47',
			cc => 'Print "error"',
			cg => '0.79',
			ga => '0.08',
			gt => '-0.52',
			gc => '-1.11',
			gg => 'Print "error"',
		);
	
	
	}elsif($Aligned eq "c"){
		my %MX_GAP_G = (
			aa => 'Print "error"',
			at => '0.43',
			ac => '0.03',
			ag => '0.75',
			ta => '-0.12',
			tt => 'Print "error"',
			tc => '-0.32',
			tg => '0.4',
			ca => '0.62',
			ct => '0.79',
			cc => 'Print "error"',
			cg => '0.7',
			ga => '-0.47',
			gt => '0.11',
			gc => '-0.11',
			gg => 'Print "error"',

		);
		$gvalue=$MX_GAP_G{$mismatch};
	}elsif($Aligned eq "t"){
		my %MX_GAP_G = (
			aa => 'Print "error"',
			at => '0.69',
			ac => '0.42',
			ag => '0.92',
			ta => '0.68',
			tt => 'Print "error"',
			tc => '0.34',
			tg => '0.75',
			ca => '0.97',
			ct => '1.33',
			cc => 'Print "error"',
			cg => '1.05',
			ga => '0.43',
			gt => '0.74',
			gc => '0.44',
			gg => 'Print "error"',

		
		);	
		$gvalue=$MX_GAP_G{$mismatch};
	}elsif($Aligned eq "a"){
		my %MX_GAP_G = (
			aa => 'Print "error"',
			at => '0.61',
			ac => '0.14',
			ag => '0.88',
			ta => '0.69',
			tt => 'Print "error"',
			tc => '0.07',
			tg => '0.73',
			ca => '0.64',
			ct => '0.77',
			cc => 'Print "error"',
			cg => '1.33',
			ga => '0.71',
			gt => '0.02',
			gc => '-0.13',
			gg => 'Print "error"',
		);
		$gvalue=$MX_GAP_G{$mismatch};
	}
	chomp $gvalue;
	return $gvalue;
	#print $gvalue."\n";
}
sub Double_Mismatch_Calculate{
	my $Query = shift;
	#print $Query."query\n";
	my $Subject =shift;
	#print $Subject."Subject\n";
	#This is the final return value variable
	my $deltag=Mismatch_calculate(substr($Query,0,2),substr($Subject,0,2));
	my $Rquery=substr($Query,2,2);
	my $RSubject=substr($Subject,2,2);
	$Rquery=reverse $Rquery;
	$RSubject=reverse $RSubject;
	$deltag=$deltag+Mismatch_calculate($RSubject,$Rquery);
	return $deltag;
}
sub Match_calculate{
	my $Query = shift;
#	print $Query."query\n";
	my $Subject =shift;
#	print $Subject."Subject\n";
	my %NO_GAP_G = (
	###the thermodynamic paramaters for an exact match. taken from: 
		at_at => '-0.88',
		ac_ac => '-1.44',
		ag_ag => '-1.28',
		ta_ta => '-0.58',
		tc_tc => '-1.3',
		tg_tg => '-1.45',
		ca_ca => '-1.45',
		ct_ct => '-1.28',
		cg_cg => '-2.17',
		ga_ga => '-1.3',
		gc_gc => '-2.24',
		gt_gt => '-1.44',	
		aa_aa => '-1',
		tt_tt => '-1',
		cc_cc => '-1.84',
		gg_gg => '-1.84',
	###the thermodynamic constants for mismatches. taken from: 
	);


	my $dinucleotide=$Query.$Subject."_".$Query.$Subject;
#	print $dinucleotide;
	return $NO_GAP_G{$dinucleotide};
##Key Variables
###gibbs free energy

	#### Enthalpy (probably won't use)
	my %NO_gaP_H = (
	###the thermodynamic paramaters for an exact match. taken from: 
		at_at => '-7.2',
		ac_ac => '-8.4',
		ag_ag => '-7.8',
		ta_ta => '-7.2',
		tc_tc => '-8.2',
		tg_tg => '-8.5',
		ca_ca => '-8.5',
		ct_ct => '-7.8',
		cg_cg => '-10.6',
		ga_ga => '-8.2',
		gc_gc => '-9.8',
		gt_gt => '-8.4',	
		aa_aa => '-7.6',
		tt_tt => '-7.6',
		cc_cc => '-8',
		gg_gg => '-8',
#	##the thermodynamic constants for mismatches. taken from: 
	);

	###Entropy (probably won't use)
	my %NO_GAP_S = (

	###the thermodynamic paramaters for an exact match. taken from: 
		at_at => '-20.4',
		ac_ac => '-22.4',
		ag_ag => '-21',
		ta_ta => '-21.3',
		tc_tc => '-22.2',
		tg_tg => '-22.7',
		ca_ca => '-22.7',
		ct_ct => '-21',
		cg_cg => '-27.2',
		ga_ga => '-22.2',
		gc_gc => '-24.4',
		gt_gt => '-22.4',	
		aa_aa => '-21.3',
		tt_tt => '-21.3',
		cc_cc => '-19',
		gg_gg => '-19',
	###the thermodynamic constants for mismatches. taken from: 
	);

}


#####Insert the query string and the subject string with one extra nucleotide in either direction.
sub Gap_Calc{
	my $Query = shift;
	#print $Query."query\n";
	my $Subject =shift;
	#print $Subject."Subject\n";
	#This is the final return value variable
	my $deltag=0;
	
	my $g1 = () = $Query =~ /-/g; 
	my $g2 = () = $Subject =~ /-/g; 
	my $AsymetryCalc=abs($g1-$g2)*0.3;

	my $n1 = () = $Query =~ /N/g; 
	my $n2 = () = $Subject =~ /N/g; 
	my $Ntotal=$n1+$n2;

	my $r1 = () = $Query =~ /R/g; 
	my $r2 = () = $Subject =~ /R/g; 
	my $Rtotal=$r1+$r2;
	
		#print length $Query;
		#print "lengthquer\n";
		#print $Ntotal."Ntotal\n";
		#print $g1."g1\n";
		#print $g2."g2\n";

	
	#Single bp bulge
	if (($g1==1 or $g2==1) and (length $Query ==3)and $Ntotal==0 and $Rtotal==0){
		$deltag=&Single_Bulge_Calc($Query,$Subject);
##		#print "SINGLE BASE PAIR BULGE\n\n\n\n\n\n";
		
		#Double bp bulge TAKEN FROM Annu. Rev. Biophys. Biomol. Struct. 2004. 33:415–40
	}elsif (($g1==2 or $g2==2) and (length($Query) ==4) and $Ntotal==0 and $Rtotal==0){
		#print "\n\nDouble bp bulge CALCULATOR";
		if((substr($Query,0,1) eq A or substr($Query,0,1) eq T) and (substr($Query,3,1) eq A or substr($Query,3,1) eq T)){
			$deltag=3.4;
		}else{
			$deltag=2.9;
		}
	}elsif (($g1==0 and $g2==0) and (length($Query) ==4) and $Ntotal==0 and $Rtotal==0){###CHECKED AND GOOD
		#### Double mismatch calculator
	$deltag=Double_Mismatch_Calculate($Query,$Subject);
	#print "Double mismatch CALCULATOR\n";

	}else{
		#Basic internal loop calculation. ASSUMING THAT THERE ISN"t AN OFFSET AND A MATCH. MAYBE ADD THIS LATER. actually it should be in the alignment unlesss its an inversion
		#TAKEN FROM Annu. Rev. Biophys. Biomol. Struct. 2004. 33:415–40
		# Remember each mismatched position counts as position.

		my $length = length($Query)+length($Subject)-$g1-$g2-4-$Ntotal;
			#print $length."length\n";

		if ($length <= 3){$length=3;}
		#print $length."length\n";
		#print log($length/30)."\n";
		#print log($length/30)."\n";
		#print $AsymetryCalc."AsymetryCalc\n";
	
		$deltag=$AsymetryCalc+6.6+2.44*0.00198788*310.15*(log($length/30));
	
	}	

	return $deltag;
}


### This one takes a three letter input with a match for positions 1 and 3. The middle position in one sequence must be a gap (-)
sub Single_Bulge_Calc{####CHECKED
	#print "\n\nBULGE CALCULATOR\n";
	my $Query = shift;
	my $Subject =shift;
	my $deltaGbulge;
	my %Bulge_G = (
	#DATA TAKEN FROM Biochemistry 2004, 43, 7143-7150
		CAA => '-0.11',
		CAT => '1.06',
		CAC => '1.73',
		CAG => '0.83',
		CTA => '0.28',
		CTT => '0.83',
		CTC => '0.57',
		CTG => '0.79',
		CCA => '1.52',
		CCT => '0.81',
		CCC => '0.51',
		CCG => '1.19',
		CGA => '0.57',
		CGT => '1',
		CGC => '1.29',
		CGG => '-0.73',
		GAA => '0.96',
		GAT => '1.19',
		GAC => '0.98',
		GAG => '0.52',
		GTA => '1.03',
		GTT => '1.99',
		GTC => '0.29',
		GTG => '0.87',
		GCA => '1.73',
		GCT => '1.42',
		GCC => '0.6',
		GCG => '2.02',
		GGA => '1.5',
		GGT => '1.75',
		GGC => '0.48',
		GGG => '0.98',
		AAA => '1.66',
		AAT => '0.35',
		AAC => '1.95',
		AAG => '1.24',
		ATA => '2.89',
		ATT => '-0.39',
		ATC => '1.17',
		ATG => '1.57',
		ACA => '2.51',
		ACT => '1.3',
		ACC => '-1.05',
		ACG => '1.07',
		AGA => '3.33',
		AGT => '1.98',
		AGC => '1.33',
		AGG => '1.29',
		TAA => '3.46',
		TAT => '3.69',
		TAC => '1.67',
		TAG => '2.36',
		TTA => '2.32',
		TTT => '0.96',
		TTC => '1.83',
		TTG => '1.38',
		TCA => '2.74',
		TCT => '2.58',
		TCC => '0.16',
		TCG => '1.32',
		TGA => '2.13',
		TGT => '3.42',
		TGC => '2.3',
		TGG => '0.33',
	);
	if ($Query =~ "N" or $Subject =~ "N"){
		$deltaGbulge=0;
	}elsif ($Query =~ "-"){

		$Subject=~tr/atcg/ATCG/;

		$deltaGbulge=$Bulge_G{$Subject};
	}else{
		$Query=~tr/atcg/ATCG/;
		$deltaGbulge=$Bulge_G{$Query};
	}
	return $deltaGbulge;
}






