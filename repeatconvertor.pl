#!/usr/bin/env perl
$file_in = $ARGV[0];
$file_in2 = $ARGV[1];
$file_out = $ARGV[2];
##this value is whether any print out is displayed
$display = $ARGV[3];
my $Disp_Out = 1;

if ($#ARGV == 2){
	$Disp_Out  = 0;
}

if ($Disp_Out  == 1){
	print "This program takes the output from Repeat convertor for a multifasta file and specifies which positions in the corresponding alignment contain gaps\n";
}

###Open files
open my $handle, '<', $file_in;
chomp(my @infile = <$handle>);
close $handle;

open my $handle, '<', $file_in2;
chomp(my @repeatlist = <$handle>);
close $handle;



###Remove headers

my @outlist;
my @sequences;
foreach my $line (@infile){
	unless ($line =~ /^>/){
		push @sequences, $line;
	}
}

my @repeats;
foreach my $line (@repeatlist){
	unless ($line =~ /^>/){
		push @repeats, $line;
	}
}




my $seqnumber;
open OF, ">", "$file_out" or die $!;

while ($seqnumber<=$#sequences){
	$sequences[$seqnumber]=~ tr/atcg/ATCG/;
	$repeats[$seqnumber]=~ tr/atcg/ATCG/;
	@nucleotides=split("",$sequences[$seqnumber]);
	@repnucs=split("",$repeats[$seqnumber]);
	my $counter =0;
	my $repadder= 0;
	while ($counter<=$#nucleotides){
		if ($nucleotides[$counter] eq "-"){$repadder++;
			# print "$repadder\nHOLY SHITBALLS\n";
		 }else{
		#	 print $nucleotides[$counter]."_".$repnucs[$counter-$repadder]."\n";
			 if($repnucs[$counter-$repadder] eq "N" and $nucleotides[$counter] ne "N"){
			 	print OF $counter.",";
			 }
		}
		


		$counter ++;
	}
	
	$seqnumber++;
	print OF "\n";
}
close OF;





