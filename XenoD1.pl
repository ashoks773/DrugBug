#!/usr/bin/perl

use strict;

#----------------------------------------XXXXXXXXXXXXXX------------------------------------------------
# This is the pipeline for the prediction of EC Class at very first level for a query molecule#
# Date: 20-Nov-2015
# Written by: Ashok K. Sharma
#----------------------------------------XXXXXXXXXXXXXX------------------------------------------------

my $inputfile = $ARGV[0]; # input molecule file (.mol) or Cid
my $model = $ARGV[1]; # input will be either a fingerprint/descripter or Hybrid
my $sample = $ARGV[2]; # input will be either normal or upsampling
my $PP = $ARGV[3]; # Additional option Of prediction probabilities
my $Thr; my $threshold; 

if ($PP == "blank") { $Thr = 0.5;}
else { $Thr = $PP; }
#print $Thr; <stdin>;

system("java -jar PaDEL-Descriptor.jar removesalt -standardizetautomers -tautomerlist tautomerlist_SMIRKS.txt -convert3d -fingerprints -maxruntime -1 -log -usefilenameasmolname maxcpdperfile 0 -dir $inputfile -file Hybrid_fingerprints.out");

#------------------------------------------------------------ Hybrid------------------------------------------------------------------------------------------------#

if ($model =~ /Hybrid/)
{	
	  system ("cut -d\",\" -f2-10209 Hybrid_fingerprints.out > Hybrid_fingerprints.Out");
	  system ("cut -d\",\" -f17,18,87,92,122,128,159,193,200,212,221,239,243,254,258,266,288,315,356,359,394,412,464,477,495,508,512,514,525,529,549,562,598,633,658,729,744,801,802,832,837,853,884,929,937,949,974,997,1002,1005,1007 Hybrid_fingerprints.Out > 1");
	  system ("cut -d\",\" -f1032,1035,1042,1045,1048,1057,1058,1060 Hybrid_fingerprints.Out > 2");
	  system ("cut -d\",\" -f1133,1195,1196,1253,1259,1269,1331,1345,1367,1389,1409,1410,1412,1431,1436,1465,1801,1825,1826,1828,1834,1855,1937,1989,2023,2085,2108 Hybrid_fingerprints.Out > 3");
	  system ("cut -d\",\" -f2140,2153,2161,2170,2176,2177,2180,2184,2211,2217,2225,2233,2237,2250,2252,2253,2263,2264,2267,2269,2273,2282,2285  Hybrid_fingerprints.Out > 4");
	  system ("cut -d\",\" -f2307,2309,2314,2315,2475,2480,2549,2617,2627,2659,2674,2699,2726,2730,2734,2745,2756,2822,2908,2953,2956,2958,2980,3112 Hybrid_fingerprints.Out > 5");
	  system ("cut -d\",\" -f3190,3192,3221,3222,3223,3249,3262,3317,3378,3449,3455,3456,3470 Hybrid_fingerprints.Out > 6");
	  system ("cut -d\",\" -f3489,3497,3670,4158,4887,5047,6028,6045,6601,6641,6749,6863,6869,7135,7149,7180,7210,7225,7270,7287,7301,7302,7363,7438,7496,7603,7674,7718,7776,8176,8292,8298,8302,8307 Hybrid_fingerprints.Out > 7");
	  system ("cut -d\",\" -f8343,8433,8443,8677,8679,8744,8911,8989 Hybrid_fingerprints.Out > 8");
	  system ("cut -d\",\" -f9124,9126,9137,9168,9169,9170,9209,9264,9325,9395,9396,9402,9403,9416,9417,9423 Hybrid_fingerprints.Out > 9");
	  system ("cut -d\",\" -f9429,9430,9431,9508,9520,9530,9585,9608,9667,9764,9831,9842,9910,9998,10076 Hybrid_fingerprints.Out > 10");

	  system ("paste -d\",\" 1 2 3 4 5 6 7 8 9 10 > Hybrid_fingerprints.Out.Final && rm 1 2 3 4 5 6 7 8 9 10 Hybrid_fingerprints.Out");

	  if ($sample =~ /Upsampling/) { system ("sh RscriptHybridUp.sh"); }
	  if ($sample =~ /Normal/) { system ("sh RscriptHybrid.sh"); }
	  system ("sed -i '1d' Hybrid_Pred");

	  my $Rfpred = "Hybrid_Pred"; # prediction file with probabilities obtained from the random forest commnad.
	  $threshold = $Thr;       # $ARGV[1]; # cutoff thresold for the prediction 
	  open (PR, $Rfpred) || die "Cant open the file:\n";
	  open (PRR, ">$Rfpred.response") || die "Cant open the file:\n";

	  while (chomp (my $line = <PR>))
	  {
	  	  # "1"	0.832	0.084	0.078	0.006	0	0
		  
		  my @arr = split ("\t", $line);	  
		  my %hash; $hash{$arr[1]}="EC1"; $hash{$arr[2]}="EC2"; $hash{$arr[3]}="EC3"; $hash{$arr[4]}="EC4"; $hash{$arr[5]}="EC5"; $hash{$arr[6]}="EC6";
		  my $x = 0; foreach my $y (@arr){ $x = $y if $y > $x; }

		  if ($arr[1] > $threshold) { print PRR $arr[0]."\tEC1\n"; }
		  elsif ($arr[2] > $threshold) { print PRR $arr[0]."\tEC2\n"; }
		  elsif ($arr[3] > $threshold) { print PRR $arr[0]."\tEC3\n"; }
		  elsif ($arr[4] > $threshold) { print PRR $arr[0]."\tEC4\n"; }
		  elsif ($arr[5] > $threshold) { print PRR $arr[0]."\tEC5\n"; }
		  elsif ($arr[6] > $threshold) { print PRR $arr[0]."\tEC6\n"; }
		  else { print PRR $arr[0]."\tNP (Probability below thresold)\t$hash{$x}\n"; } # We have to decide wheteher this class should included or not 
	  }
close PR; close PRR;
}
#------------------------------------------------ Fingerprinter----------------------------------------------------------------------------------------------------#
if ($model =~ /Fingerprinter/)
{
	  system ("cut -d\",\" -f2-1025 Hybrid_fingerprints.out > Fingerprinter.Out.Final");
	  system ("rm Fingerprinter.out");

	  if ($sample =~ /Upsampling/){ system ("sh RscriptFingerprinterUp.sh"); }
	  if ($sample =~ /Normal/) { system ("sh RscriptFingerprinter.sh"); }
	  system ("sed -i '1d' Fingerprinter_Pred");

	  my $Rfpred = "Fingerprinter_Pred"; # prediction file with probabilities obtained from the random forest commnad.
	  $threshold = $Thr;    # $ARGV[1]; # cutoff thresold for the prediction 
	  open (PR, $Rfpred) || die "Cant open the file:\n";
	  open (PRR, ">$Rfpred.response") || die "Cant open the file:\n";

	  while (chomp (my $line = <PR>))
	  {
	  # "1"	0.832	0.084	0.078	0.006	0	0
	  # "2"	0.148	0.586	0.16	0.086	0.014	0.006
		  
		  my @arr = split ("\t", $line);
		  my %hash; $hash{$arr[1]}="EC1"; $hash{$arr[2]}="EC2"; $hash{$arr[3]}="EC3"; $hash{$arr[4]}="EC4"; $hash{$arr[5]}="EC5"; $hash{$arr[6]}="EC6";
		  my $x = 0; foreach my $y (@arr){ $x = $y if $y > $x; }

		  if ($arr[1] > $threshold) { print PRR $arr[0]."\tEC1\n"; }
		  elsif ($arr[2] > $threshold) { print PRR $arr[0]."\tEC2\n"; }
		  elsif ($arr[3] > $threshold) { print PRR $arr[0]."\tEC3\n"; }
		  elsif ($arr[4] > $threshold) { print PRR $arr[0]."\tEC4\n"; }
		  elsif ($arr[5] > $threshold) { print PRR $arr[0]."\tEC5\n"; }
		  elsif ($arr[6] > $threshold) { print PRR $arr[0]."\tEC6\n"; }
		  else { print PRR $arr[0]."\tNP (Probability below thresold)\t$hash{$x}\n"; } # We have to decide wheteher this class should included or not 
	  }
close PR; close PRR;
}
#----------------------------------------- Pubchem-----------------------------------------------------------------------------------------------------------------#
if ($model =~ /Pubchem/)
{
	  system ("cut -d\",\" -f2295-3176 Hybrid_fingerprints.out > Pubchem.Out.Final");
	  system ("rm Pubchem.out");

	  if ($sample =~ /Upsampling/) { system ("sh RscriptPubchemUp.sh"); }
	  if ($sample =~ /Normal/) {  system ("sh RscriptPubchem.sh"); }
	  system ("sed -i '1d' Pubchem_Pred");

	  my $Rfpred = "Pubchem_Pred"; # prediction file with probabilities obtained from the random forest commnad.
	  $threshold = $Thr;        # $ARGV[1]; # cutoff thresold for the prediction 
	  open (PR, $Rfpred) || die "Cant open the file:\n";
	  open (PRR, ">$Rfpred.response") || die "Cant open the file:\n";

	  while (chomp (my $line = <PR>))
	  {
	  # "1"	0.832	0.084	0.078	0.006	0	0
	  # "2"	0.148	0.586	0.16	0.086	0.014	0.006
		  
		  my @arr = split ("\t", $line);
		  my %hash; $hash{$arr[1]}="EC1"; $hash{$arr[2]}="EC2"; $hash{$arr[3]}="EC3"; $hash{$arr[4]}="EC4"; $hash{$arr[5]}="EC5"; $hash{$arr[6]}="EC6";
		  my $x = 0; foreach my $y (@arr){ $x = $y if $y > $x; }

		  if ($arr[1] > $threshold) { print PRR $arr[0]."\tEC1\n"; }
		  elsif ($arr[2] > $threshold) { print PRR $arr[0]."\tEC2\n"; }
		  elsif ($arr[3] > $threshold) { print PRR $arr[0]."\tEC3\n"; }
		  elsif ($arr[4] > $threshold) { print PRR $arr[0]."\tEC4\n"; }
		  elsif ($arr[5] > $threshold) { print PRR $arr[0]."\tEC5\n"; }
		  elsif ($arr[6] > $threshold) { print PRR $arr[0]."\tEC6\n"; }
		  else { print PRR $arr[0]."\tNP (Probability below thresold)\t$hash{$x}\n"; } # We have to decide wheteher this class should included or not 
	  }
close PR; close PRR;
}


